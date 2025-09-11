// Elements
const fileInput = document.getElementById('ecgFile');
const analyzeBtn = document.getElementById('analyzeBtn');
const resultText = document.getElementById('resultText');
const canvas = document.getElementById('ecgCanvas');
const ctx = canvas.getContext('2d');

// Setup canvas size
canvas.width = window.innerWidth * 0.8;
canvas.height = 300;

// Dummy ECG generator
function generateDummyECG() {
  const points = [];
  for (let x = 0; x < canvas.width; x++) {
    const y = 150 + Math.sin(x * 0.05) * 50 + Math.random() * 10;
    points.push(y);
  }
  return points;
}

// Draw ECG waveform
function drawECG(points) {
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.strokeStyle = "#00ffcc";
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.moveTo(0, points[0]);
  points.forEach((y, x) => ctx.lineTo(x, y));
  ctx.stroke();
}

// Fake analyzer
function analyzeECG() {
  const outcomes = [
    "✅ Normal Sinus Rhythm",
    "⚠️ Possible Atrial Fibrillation",
    "⚠️ Bradycardia Detected",
    "⚠️ Tachycardia Detected",
    "❌ ECG Abnormal"
  ];
  const result = outcomes[Math.floor(Math.random() * outcomes.length)];
  resultText.textContent = result;

  // Show ECG waveform
  drawECG(generateDummyECG());
}

// Click listener
analyzeBtn.addEventListener('click', () => {
  if (fileInput.files.length > 0) {
    analyzeECG();
  } else {
    resultText.textContent = "Please upload an ECG file first!";
  }
});
