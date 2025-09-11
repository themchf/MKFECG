/**
 * MKF-ECG — pure JS high-precision signal pipeline
 *  - bandpass FIR (sinc window)
 *  - derivative -> squaring -> moving-window integration (Pan-Tompkins style)
 *  - adaptive thresholding for QRS detection
 *  - refine R-peak location on filtered signal
 *  - compute RR intervals and HRV metrics (SDNN, RMSSD, pNN50)
 *  - rule-based detection: brady/tachy, AF-likelihood, PVCs, long pause
 *
 * Usage: upload CSV (single column). Default fs=250 Hz.
 */

/* ---------- helpers ---------- */

// convolution (direct)
function convolve(signal, kernel) {
  const n = signal.length, m = kernel.length;
  const out = new Float32Array(n);
  const kHalf = Math.floor(m/2);
  for (let i=0;i<n;i++){
    let s = 0;
    for (let j=0;j<m;j++){
      const idx = i - j + kHalf;
      if (idx >=0 && idx < n) s += signal[idx] * kernel[j];
    }
    out[i] = s;
  }
  return out;
}

// design FIR bandpass via windowed sinc (Hamming)
function firBandpass(lowcut, highcut, fs, N=101) {
  if (N%2===0) N++; // odd
  const fc1 = lowcut / fs;
  const fc2 = highcut / fs;
  const h = new Float32Array(N);
  const M = (N-1)/2;
  for (let n=0;n<N;n++){
    const k = n - M;
    if (k === 0) {
      h[n] = 2*(fc2 - fc1);
    } else {
      h[n] = (Math.sin(2*Math.PI*fc2*k) - Math.sin(2*Math.PI*fc1*k)) / (Math.PI * k);
    }
    // Hamming window
    h[n] *= 0.54 - 0.46 * Math.cos(2*Math.PI*n/(N-1));
  }
  return h;
}

// moving average (boxcar)
function movingAverage(arr, win) {
  const n=arr.length; const out=new Float32Array(n);
  let sum=0;
  for (let i=0;i<n;i++){
    sum += arr[i];
    if (i>=win) sum -= arr[i-win];
    out[i] = sum/Math.min(i+1, win);
  }
  return out;
}

// standard deviation
function std(arr) {
  const n=arr.length;
  if(n===0) return 0;
  const mean = arr.reduce((a,b)=>a+b,0)/n;
  const v = arr.reduce((a,b)=>a+((b-mean)*(b-mean)),0)/n;
  return Math.sqrt(v);
}

/* ---------- Pan-Tompkins style pipeline ---------- */

function preprocessAndDetect(raw, fs) {
  // 1) detrend baseline by subtracting long moving average (0.6s)
  const baselineWin = Math.round(0.6 * fs);
  const baseline = movingAverage(raw, baselineWin);
  const detr = new Float32Array(raw.length);
  for (let i=0;i<raw.length;i++) detr[i] = raw[i] - baseline[i];

  // 2) bandpass for QRS enhancement (5–15 Hz typical)
  const bp = firBandpass(5, 15, fs, Math.max(101, Math.round(fs*0.2)|1)); // N ~ 0.2*samples/sec
  const filtered = convolve(detr, bp);

  // 3) derivative (5-point)
  const deriv = new Float32Array(filtered.length);
  for (let i=2;i<filtered.length-2;i++){
    deriv[i] = (2*filtered[i+2] + filtered[i+1] - filtered[i-1] - 2*filtered[i-2]) / 8.0;
  }

  // 4) squaring
  const squared = new Float32Array(deriv.length);
  for (let i=0;i<deriv.length;i++) squared[i] = deriv[i]*deriv[i];

  // 5) moving window integration (window ~ 150 ms)
  const intWin = Math.max(1, Math.round(0.15 * fs));
  const integrated = movingAverage(squared, intWin);

  // 6) adaptive thresholding & locate QRS on integrated signal
  const peaks = detectPeaksAdaptive(integrated, filtered, fs);

  // 7) compute RR intervals and metrics
  const rrsec = [];
  for (let i=1;i<peaks.length;i++){
    rrsec.push((peaks[i] - peaks[i-1]) / fs);
  }
  const meanRR = rrsec.length? (rrsec.reduce((a,b)=>a+b,0)/rrsec.length):0;
  const sdnn = rrsec.length? std(rrsec):0;
  const diffs = [];
  for (let i=1;i<rrsec.length;i++) diffs.push(Math.abs(rrsec[i]-rrsec[i-1]));
  const rmssd = Math.sqrt(diffs.length? (diffs.reduce((a,b)=>a+b*b,0)/diffs.length):0);
  const pnn50 = diffs.length? (diffs.filter(d=>d>0.05).length / diffs.length):0; // >50ms -> 0.05s

  return {
    filtered, integrated, peaks,
    rrsec, meanRR, sdnn, rmssd, pnn50
  };
}

/* Adaptive detection: find peaks in integrated signal, refine on filtered */
function detectPeaksAdaptive(integrated, filtered, fs) {
  // Find local maxima in integrated
  const N = integrated.length;
  const localPeaks = [];
  for (let i=2;i<N-2;i++){
    if (integrated[i] > integrated[i-1] && integrated[i] > integrated[i+1]) {
      localPeaks.push(i);
    }
  }
  if (localPeaks.length === 0) return [];

  // Adaptive threshold based on signal & noise estimation
  // Use initial 2 seconds to seed thresholds if possible
  const seedLen = Math.min(N, Math.round(fs*2));
  const seedPeaks = localPeaks.filter(p => p < seedLen);
  let SPKI = seedPeaks.length? (seedPeaks.reduce((a,b)=>a+integrated[b],0)/seedPeaks.length) : Math.max(...integrated)*0.3;
  let NPKI = 0.02 * SPKI;
  let threshI1 = NPKI + 0.25*(SPKI - NPKI);

  const qrsCandidates = [];
  let lastQRS = -9999;
  const refractory = Math.round(0.25 * fs); // 250ms refractory

  for (const idx of localPeaks) {
    const val = integrated[idx];
    if (val > threshI1 && (idx - lastQRS) > refractory) {
      // refine R position search in filtered signal around idx ± 50 ms
      const searchHalf = Math.round(0.05*fs);
      const start = Math.max(0, idx-searchHalf);
      const end = Math.min(filtered.length-1, idx+searchHalf);
      let rIdx = start;
      let rVal = filtered[start];
      for (let j=start+1;j<=end;j++){
        if (filtered[j] > rVal) { rVal = filtered[j]; rIdx = j; }
      }
      qrsCandidates.push(rIdx);
      lastQRS = idx;
      // update signal level
      SPKI = 0.125*val + 0.875*SPKI;
    } else {
      // update noise level
      NPKI = 0.125*val + 0.875*NPKI;
    }
    threshI1 = NPKI + 0.25*(SPKI - NPKI);
  }

  // Remove duplicates/close detections: keep only one per 200 ms
  const finalPeaks=[];
  for (let p of qrsCandidates){
    if (finalPeaks.length===0 || (p - finalPeaks[finalPeaks.length-1]) > Math.round(0.2*fs)) finalPeaks.push(p);
  }

  return finalPeaks;
}

/* ---------- Rule engine / medical heuristics ---------- */

function interpretFindings(rrsec, meanRR, sdnn, rmssd, pnn50, peaks, fs) {
  const out = [];
  const bpm = meanRR? Math.round(60/meanRR) : 0;

  // Heart rate categories
  if (bpm === 0) out.push("No reliable heartbeats detected");
  else {
    if (bpm < 50) out.push(`Bradycardia (HR=${bpm} bpm) — clinically significant`);
    else if (bpm < 60) out.push(`Mild bradycardia (HR=${bpm} bpm)`);
    else if (bpm <= 100) out.push(`Normal heart rate (HR=${bpm} bpm)`);
    else if (bpm <= 150) out.push(`Tachycardia (HR=${bpm} bpm)`);
    else out.push(`High rate (HR=${bpm} bpm) — consider urgent evaluation`);
  }

  // AF detection heuristic (rule-based):
  // Use combination of SDNN, RMSSD and pNN50 thresholds (seconds)
  const afRules = [
    sdnn > 0.12,        // SDNN > 120 ms
    rmssd > 0.1,        // RMSSD > 100 ms
    pnn50 > 0.20        // pNN50 > 20%
  ];
  const afCount = afRules.filter(x=>x).length;
  if (afCount >= 2) out.push("Rhythm: Likely AFib (irregularly irregular)");
  else if (afCount === 1) out.push("Rhythm: Possible AFib — correlate clinically");
  else out.push("Rhythm: No strong AFib signature detected");

  // PVC detection: look for premature beat followed by compensatory pause
  const pvcIndices = [];
  if (rrsec.length >= 3) {
    const mean = meanRR;
    for (let i=1;i<rrsec.length-1;i++){
      const prev=rrsec[i-1], curr=rrsec[i], next=rrsec[i+1];
      // premature: curr < 0.8*prev and next > 1.15*prev (compensatory pause)
      if (curr < 0.8*prev && next > 1.15*prev) {
        pvcIndices.push(i); // mark beat index (the premature one)
      }
    }
  }
  if (pvcIndices.length>0) out.push(`PVCs detected: ${pvcIndices.length} candidate(s)`);
  // long pause
  const pauses = rrsec.filter(r=>r > 3.0);
  if (pauses.length>0) out.push(`Long pause(s) detected (>3s) : ${pauses.length}`);

  // HRV summary
  out.push(`HRV: SDNN=${(sdnn||0).toFixed(3)}s, RMSSD=${(rmssd||0).toFixed(3)}s, pNN50=${Math.round((pnn50||0)*100)}%`);

  return { text: out.join(" · "), details: { bpm, sdnn, rmssd, pnn50, pvcCount: pvcIndices.length, longPauses: pauses.length } };
}

/* ---------- plotting helpers ---------- */
let ecgChart = null, rrChart=null;
function plotECGCanvas(filtered, peaks, fs) {
  const canvas = document.getElementById('ecgChart');
  const n = filtered.length;
  const labels = new Array(n);
  for (let i=0;i<n;i++) labels[i] = (i/fs).toFixed(2);
  const peakPoints = new Array(n).fill(null);
  for (let p of peaks) peakPoints[p] = filtered[p];

  const ctx = canvas.getContext('2d');
  // use Chart.js for nicer visuals
  if (ecgChart) ecgChart.destroy();
  ecgChart = new Chart(ctx, {
    type:'line',
    data:{
      labels: labels,
      datasets:[
        { label:'Filtered ECG', data: Array.from(filtered), borderColor:'#00ffd1', borderWidth:1, pointRadius:0 },
        { label:'R peaks', data: peakPoints, borderColor:'#ff3b6b', borderWidth:0, pointRadius:4, showLine:false }
      ]
    },
    options:{
      animation:false,
      scales:{ x:{ display:true, title:{display:true,text:'Time (s)'} }, y:{ display:true, title:{display:true,text:'Amplitude'} } },
      elements:{ line:{ tension:0 } },
      plugins:{ legend:{ position:'top' } }
    }
  });
}

function plotRRHistogram(rrsec) {
  const canvas = document.getElementById('rrChart');
  const ctx = canvas.getContext('2d');
  // build histogram bins (0.2s bins)
  const bins = {};
  for (let r of rrsec) {
    const k = Math.round(r*5)/5; // 0.2s steps
    bins[k] = (bins[k]||0) + 1;
  }
  const ks = Object.keys(bins).sort((a,b)=>a-b);
  const vals = ks.map(k=>bins[k]);
  if (rrChart) rrChart.destroy();
  rrChart = new Chart(ctx, {
    type:'bar',
    data:{ labels: ks, datasets:[{ label:'RR interval (s) counts', data: vals, backgroundColor:'#00aaff'}] },
    options:{ animation:false, scales:{ x:{ title:{display:true,text:'RR (s)'} }, y:{ title:{display:true,text:'Count'} } }, plugins:{ legend:{ display:false } } }
  });
}

/* ---------- UI / wiring ---------- */

document.getElementById('analyzeBtn').addEventListener('click', ()=> {
  const f = document.getElementById('ecgFile').files[0];
  if(!f){ alert('Upload a single-column CSV containing ECG samples (voltage)'); return; }
  const fs = Number(document.getElementById('fs').value) || 250;
  const winLen = Number(document.getElementById('winLen').value) || 30;

  const reader = new FileReader();
  reader.onload = (e)=>{
    // read numeric values (strip empty lines)
    const lines = e.target.result.split(/\r?\n/).map(s=>s.trim()).filter(s=>s!=="");
    // if CSV has two columns (time,value) try to detect column with numbers length>1
    const vals = [];
    for (let L of lines){
      const cols = L.split(',');
      // pick first numeric column
      let pick = null;
      for (let c of cols){
        const v = parseFloat(c);
        if (!isNaN(v)){ pick = v; break; }
      }
      if (pick!==null) vals.push(pick);
    }
    if (vals.length < 50) { alert('ECG too short (need at least ~50 samples)'); return; }

    // optionally crop to window length
    const maxSamples = Math.round(winLen * fs);
    const raw = Float32Array.from(vals.slice(0, maxSamples));

    // run pipeline
    const result = preprocessAndDetect(raw, fs);

    // compute metrics
    const { rrsec, meanRR, sdnn, rmssd, pnn50, peaks, filtered } = result;
    // update UI
    document.getElementById('beats').innerText = peaks.length || 0;
    document.getElementById('hr').innerText = meanRR ? Math.round(60/meanRR) : '--';
    document.getElementById('sdnn').innerText = sdnn ? sdnn.toFixed(3) : '--';
    document.getElementById('rmssd').innerText = rmssd ? rmssd.toFixed(3) : '--';
    document.getElementById('pnn50').innerText = pnn50? (Math.round(pnn50*100)+'%') : '--';
    document.getElementById('meanRR').innerText = meanRR? meanRR.toFixed(3) : '--';

    const interp = interpretFindings(rrsec, meanRR, sdnn, rmssd, pnn50, peaks, fs);
    document.getElementById('summary').innerText = interp.text;

    plotECGCanvas(filtered, peaks, fs);
    plotRRHistogram(rrsec);
  };
  reader.readAsText(f);
});
