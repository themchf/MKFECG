const analyzeBtn = document.getElementById('analyzeBtn');
const ecgFile = document.getElementById('ecgFile');
const output = document.getElementById('output');
const ecgChart = document.getElementById('ecgChart').getContext('2d');

// CSV parser
function parseCSV(file, callback) {
    const reader = new FileReader();
    reader.onload = e => {
        const lines = e.target.result.trim().split('\n');
        const data = lines.map(line => parseFloat(line.trim())).filter(x => !isNaN(x));
        callback(data);
    };
    reader.readAsText(file);
}

// Normalize data to 0-1
function normalize(data) {
    const min = Math.min(...data);
    const max = Math.max(...data);
    return data.map(x => (x - min) / (max - min));
}

// Peak detection with dynamic threshold
function detectPeaks(data) {
    const peaks = [];
    const threshold = 0.6 * Math.max(...data);
    for (let i = 1; i < data.length - 1; i++) {
        if (data[i] > threshold && data[i] > data[i-1] && data[i] > data[i+1]) {
            peaks.push(i);
        }
    }
    return peaks;
}

// Heart rate calculation
function calculateHR(peaks, samplingRate = 250) {
    if (peaks.length < 2) return 0;
    const rrIntervals = [];
    for (let i = 1; i < peaks.length; i++) {
        rrIntervals.push((peaks[i] - peaks[i-1]) / samplingRate);
    }
    const avgRR = rrIntervals.reduce((a,b)=>a+b,0)/rrIntervals.length;
    return 60 / avgRR;
}

// Arrhythmia detection based on RR interval variability
function analyzeECG(data) {
    const peaks = detectPeaks(data);
    const bpm = calculateHR(peaks);
    let diagnosis = 'Normal';

    if (bpm === 0) diagnosis = 'No detectable signal';
    else if (bpm < 60) diagnosis = 'Bradycardia';
    else if (bpm > 100) diagnosis = 'Tachycardia';

    // Detect AFib / irregular rhythm
    const rrIntervals = [];
    for (let i = 1; i < peaks.length; i++) {
        rrIntervals.push(peaks[i] - peaks[i-1]);
    }
    const rrMean = rrIntervals.reduce((a,b)=>a+b,0)/rrIntervals.length;
    const rrStd = Math.sqrt(rrIntervals.map(x=>Math.pow(x-rrMean,2)).reduce((a,b)=>a+b,0)/rrIntervals.length);
    if (rrStd > 15) diagnosis += ' / Possible AFib';

    return { bpm: bpm.toFixed(1), diagnosis };
}

// Plot ECG using Chart.js
function plotECG(data) {
    new Chart(ecgChart, {
        type: 'line',
        data: {
            labels: data.map((_,i)=>i),
            datasets: [{
                label: 'ECG Signal',
                data: data,
                borderColor: '#ff0066',
                fill: false,
                borderWidth: 2,
                pointRadius: 0
            }]
        },
        options: {
            animation: false,
            responsive: true,
            scales: {
                x: { display: true, title: { display: true, text: 'Sample Index', color: '#00ffcc' } },
                y: { display: true, title: { display: true, text: 'Amplitude', color: '#00ffcc' } }
            },
            plugins: {
                legend: { display: false }
            }
        }
    });
}

analyzeBtn.addEventListener('click', () => {
    const file = ecgFile.files[0];
    if (!file) return alert('Please upload a CSV ECG file');

    parseCSV(file, (data) => {
        const normalized = normalize(data);
        plotECG(normalized);
        const result = analyzeECG(normalized);
        output.innerHTML = `Heart Rate: <strong>${result.bpm} BPM</strong><br>Diagnosis: <strong>${result.diagnosis}</strong>`;
    });
});
