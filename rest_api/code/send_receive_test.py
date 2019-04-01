from matchable_observable import simAlert
from requests import put, get

windowAwake = [16, 15.5, 16, 13]
alertData = simAlert(4) #number_days
alertData.generate(windowAwake, 100000) #window_decision, resolution
samples = alertData.randSample(2, 5)

put('http://localhost:5000/al_01', data=samples).json()

get('http://localhost:5000/al_01').json()
