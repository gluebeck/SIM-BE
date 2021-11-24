detectProb = function(m,scaleDetection) {
  probs = (m/scaleDetection)^1.3; probs[probs > 1]=1
  return(probs)
}
