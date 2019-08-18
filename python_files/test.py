from flask import Flask, request, jsonify
from flask_restful import Resource, Api
from flask_request_params import bind_request_params
from matchable_observable import MOLIwhiteBox, forceEstimate
import json
import pyrebase
import datetime

# Initialize the restful API app
app = Flask(__name__)
api = Api(app)
app.before_request(bind_request_params)

# Initialize some estimate models
modelMOLI = MOLIwhiteBox('trivial')
modelAnnealing = forceEstimate('Annealing')

# Initialize connection with Firebase
config = {
  "apiKey": "AIzaSyAoVwP8yx9OKWfm89MECLUiLwmnwBJvOHQ",
  "authDomain": "alertnessalpha01.firebaseapp.com",
  "databaseURL":"https://alertnessalpha01.firebaseio.com",
  "storageBucket": "alertnessalpha01.appspot.com"
}

firebase = pyrebase.initialize_app(config)
fireDB = firebase.database()

class ToDoSimple(Resource):
    def get(self, model_id):
        requested_data = request.params
        sTime = {'init': json.loads(requested_data['init']),
                 'final': json.loads(requested_data['final'])}

        __init = json.loads(requested_data['lev'])
        data = modelMOLI.predict(sTime, __init)
        return data

    def put(self, model_id):
        samples = {'time': json.loads(request.form['time']),
                   'levAl': json.loads(request.form['levAl']),
                   'windDecisions': json.loads(request.form['windDecisions'])}
        modelMOLI.fit(samples)
        return {'omega': modelMOLI.w}

class resAnnealing(Resource):
    def get(self, model_id):
        requested_data = request.params
        samples = {'initial': json.loads(requested_data['initial']),
                   'final': json.loads(requested_data['final'])}

        __initAl = json.loads(requested_data['__initAl'])

        data = modelAnnealing.simulate(samples, __initAl)

        return data

    def put(self, model_id):
        samples = {'time': json.loads(request.form['time']),
                   'levAl': json.loads(request.form['levAl']),
                   'initial': json.loads(request.form['initial']),
                   'final': json.loads(request.form['final'])}

        modelAnnealing.fit(samples)

        return {'omega': modelAnnealing.omega}

class resEstimate(Resource):
    def get(self):
        requested_data = request.params
        user_code = requested_data['user']
        # Retrive data from FireBase
        JSONsamples = fireDB.child("users").child(user_code).get().val()
        samples = {'time': json.loads(JSONsamples['time']),
                   'levAl': json.loads(JSONsamples['levAl']),
                   'initial': json.loads(JSONsamples['initial']),
                   'final': json.loads(JSONsamples['final'])}
        # Estimate the model
        modelAnnealing.fit(samples)
        # Determine the parameters
        parInfo = {'omega': modelAnnealing.omega,
                   'tau': modelAnnealing.tau,
                   'phi': modelAnnealing.phi,
                   'M': modelAnnealing.M,
                   'DC': modelAnnealing.DC,
                   'y_': modelAnnealing.y_,
                   'tau_e': modelAnnealing.tau_e}
        # Send data to Firebase
        fireDB.child("parameters").child(user_code).set(parInfo)
        # Send debug Info
        return parInfo

class fireBaseUpdate(Resource):
    def get(self):
        requested_data = request.params
        user_code = requested_data['user']
        # Retrieve last posted data from FireBase
        FireLast = fireDB.child("last").child(user_code).get().val()
        examsCond = json.loads(FireLast['exam'])
        # Determine the current time for further necessities
        currentTime = datetime.datetime.now()
        #currentTime = currentTime.replace(hour = currentTime.hour - 3)
        currentTime = currentTime - datetime.timedelta(hours=3)
        # Check the results for Classification
        if (examsCond[0]):
            # Retrieve existing data of user
            JSONsamples = fireDB.child("users").child(user_code).child("class").get().val()
            if JSONsamples != None: # if there is data already
                # Format the data into a dictionary
                FireSample = {'reference': JSONsamples['reference'],
                              'time': json.loads(JSONsamples['time']),
                              'acum_time': json.loads(JSONsamples['acum_time']),
                              'classes': json.loads(JSONsamples['classes']),
                              'before': JSONsamples['before'],
                              'lastTimeRef': JSONsamples['lastTimeRef'],
                              'missing': JSONsamples['missing'],
                              'measures': JSONsamples['measures'],
                              'Working': JSONsamples['Working'],
                              'Eating': JSONsamples['Eating'],
                              'Exercising': JSONsamples['Exercising'],
                              'Resting': JSONsamples['Resting'],
                              'Bed': JSONsamples['Bed'],
                              'Recreation': JSONsamples['Recreation'],
                              'SoftWork': JSONsamples['SoftWork'],
                              'initial': json.loads(JSONsamples['initial']),
                              'final': json.loads(JSONsamples['final'])}
                # Determine the datetime difference since reference
                form = "%d-%b-%Y (%H:%M:%S.%f)"
                datRef = datetime.datetime.strptime(FireSample['reference'], form)
                todayDiff = currentTime - datRef
                # Get the difference in hours
                days, seconds = todayDiff.days, todayDiff.seconds
                todayRef = days * 24 + seconds / 3600
                # Check the number of missing days
                # Check the number of missing days
                #if days - int(FireSample['lastTimeRef']) > 1:
                #    FireSample['missing'] = FireSample['missing'] + 1
                dist = int(days - int(FireSample['lastTimeRef']) - int(FireSample['missing']))
                if dist > 1:
                    FireSample['missing'] = FireSample['missing'] + dist - 1
                # atribute the correct day index
                index = str(days - FireSample['missing'])
                # Include the class into classes
                if FireSample['time'].get(index) == None:
                    # Determine initial vector
                    if currentTime.hour <= 6:
                        FireSample['initial'].append(todayRef)
                    else:
                        FireSample['initial'].append(todayRef - 1.5)
                    # Determine the sleeping moment
                    lastInd = str(int(index) - 1)
                    FireSample['final'].append(FireSample['time'][lastInd][-1])
                    # Append the info into the time serie
                    FireSample['classes'][index] = []
                    FireSample['acum_time'][index] = []
                    FireSample['time'][index] = [todayRef]
                else:
                    # Append the info into the time serie
                    FireSample['classes'][index].append(FireSample['before'])
                    FireSample['time'][index].append(todayRef)
                    timeGrab = todayRef - FireSample['time'][index][-2]
                    FireSample['acum_time'][index].append(timeGrab)
                    #
                    before = FireSample['before']
                    FireSample[before] = timeGrab + FireSample[before]
                # Json dumps the listsa
                FireInfo = {'reference': JSONsamples['reference'],
                            'time': json.dumps(FireSample['time']),
                            'acum_time': json.dumps(FireSample['acum_time']),
                            'classes': json.dumps(FireSample['classes']),
                            'before': FireLast['classe'],
                            'measures': int(FireSample['measures']) + 1,
                            'missing': FireSample['missing'],
                            'lastTimeRef': index,
                            'Working': FireSample['Working'],
                            'Eating': FireSample['Eating'],
                            'Exercising': FireSample['Exercising'],
                            'Resting': FireSample['Resting'],
                            'Bed': FireSample['Bed'],
                            'Recreation': FireSample['Recreation'],
                            'SoftWork': FireSample['SoftWork'],
                            'initial': json.dumps(FireSample['initial']),
                            'final': json.dumps(FireSample['final'])}
            else:
                # Initializing the first reference at 04:00:00 AM
                timeRef = currentTime.replace(hour = 4)
                timeRef = timeRef.replace(minute = 0)
                timeRef = timeRef.replace(second = 0)
                timeRef = timeRef.replace(microsecond = 0)
                # Determine the proper time reference in string
                dateInStr = timeRef.strftime("%d-%b-%Y (%H:%M:%S.%f)")
                # Determine the relative initial time
                diff = currentTime - timeRef
                # Get the difference in hours
                days, seconds = diff.days, diff.seconds
                initialRef = days * 24 + seconds / 3600
                # Create the initial time and alert dict
                __time = {"0": [initialRef]}
                aux = {"0": []}
                FireInfo = {'reference': dateInStr,
                            'time': json.dumps(__time),
                            'acum_time': json.dumps(aux),
                            'classes': json.dumps(aux),
                            'before': FireLast["classe"],
                            'lastTimeRef': 0,
                            'measures': 1, 'missing': 0,
                            'Working': 0.0, 'Eating': 0.0,
                            'Exercising': 0.0, 'Resting': 0.0,
                            'Bed': 0.0, 'Recreation': 0.0,
                            'SoftWork': 0.0,
                            'initial': json.dumps([initialRef]),
                            'final': json.dumps([])}
            # Set the data back in firebase
            fireDB.child("users").child(user_code).child("class").set(FireInfo)

        if (examsCond[1]):
            # Retrieve data samples from fireBase
            JSONsamples = fireDB.child("users").child(user_code).child("kss").get().val()
            # Transfotm sleepiness(score) into alertness
            scoreToAlert = 13*(10 - float(FireLast['score']))/9 + 1
            # Check if there is info on KSS for [user_code]
            if JSONsamples != None:
                # Format data as dictionary
                FireSample = {'time': json.loads(JSONsamples['time']),
                              'alert': json.loads(JSONsamples['alert']),
                              'reference': JSONsamples['reference'],
                              'lastTimeRef': JSONsamples['lastTimeRef'],
                              'missing': JSONsamples['missing'],
                              'measures': JSONsamples['measures'],
                              'initial': json.loads(JSONsamples['initial']),
                              'final': json.loads(JSONsamples['final'])}
                # Determine the hours expent since reference
                form = "%d-%b-%Y (%H:%M:%S.%f)"
                datRef = datetime.datetime.strptime(FireSample['reference'], form)
                todayDiff = currentTime - datRef
                # Get the difference in hours
                days, seconds = todayDiff.days, todayDiff.seconds
                todayRef = days * 24 + seconds / 3600
                # Check the number of missing days
                #if days - int(FireSample['lastTimeRef']) > 1:
                #    FireSample['missing'] = FireSample['missing'] + 1
                dist = int(days - int(FireSample['lastTimeRef']) - int(FireSample['missing']))
                if dist > 1:
                    FireSample['missing'] = FireSample['missing'] + dist - 1
                # Determine the day index to append
                index = str(days - FireSample['missing'])
                # Check if there is information on the index
                if FireSample['alert'].get(index) == None:
                    # Determine initial value vector
                    if currentTime.hour <= 6:
                        FireSample['initial'].append(todayRef)
                    else:
                        FireSample['initial'].append(todayRef - 1.5)
                    # Determine the sleeping moment
                    lastInd = str(int(index) - 1)
                    FireSample['final'].append(FireSample['time'][lastInd][-1])
                    # Include into alert
                    FireSample['alert'][index] = [scoreToAlert]
                    # Append the info into the time serie
                    FireSample['time'][index] = [todayRef]
                else:
                    # Include into alert
                    FireSample['alert'][index].append(scoreToAlert)
                    # Append the info into the time serie
                    FireSample['time'][index].append(todayRef)
                # Json dumps the lists
                FireInfo = {'time': json.dumps(FireSample['time']),
                            'alert': json.dumps(FireSample['alert']),
                            'reference': JSONsamples['reference'],
                            'measures': float(FireSample['measures']) + 1,
                            'lastTimeRef': index,
                            'missing': FireSample['missing'],
                            'initial': json.dumps(FireSample['initial']),
                            'final': json.dumps(FireSample['final'])}
            else:
                # Initializing the first reference at 04:00:00 AM
                timeRef = currentTime.replace(hour = 4)
                timeRef = timeRef.replace(minute = 0)
                timeRef = timeRef.replace(second = 0)
                timeRef = timeRef.replace(microsecond = 0)
                # Determine the proper time reference in string
                dateInStr = timeRef.strftime("%d-%b-%Y (%H:%M:%S.%f)")
                # Determine the relative initial time
                diff = currentTime - timeRef
                # Get the difference in hours
                days, seconds = diff.days, diff.seconds
                initialRef = days * 24 + seconds / 3600
                # Create the initial time and alert dict
                __alert = {"0": [scoreToAlert]}
                __time = {"0":[initialRef]}
                FireInfo = {'time': json.dumps(__time),
                            'alert': json.dumps(__alert),
                            'initial': json.dumps([initialRef]),
                            'final': json.dumps([]),
                            'measures': 1, 'missing': 0,
                            'lastTimeRef': 0,
                            'reference': dateInStr}
            fireDB.child("users").child(user_code).child("kss").set(FireInfo)

        if (examsCond[2]):
            # Retrieve data samples from fireBase
            JSONsamples = fireDB.child("users").child(user_code).child("pvt").get().val()
            # Check if there is info on PVT for [user_code]
            if JSONsamples != None:
                # Format data as dictionary
                FireSample = {'time': json.loads(JSONsamples['time']),
                              'react': json.loads(JSONsamples['react']),
                              'error': json.loads(JSONsamples['error']),
                              'initial': json.loads(JSONsamples['initial']),
                              'final': json.loads(JSONsamples['final']),
                              'reference': JSONsamples['reference'],
                              'missing': JSONsamples['missing'],
                              'lastTimeRef': JSONsamples['lastTimeRef'],
                              'measures': JSONsamples['measures']}
                # Determine the hours expent since reference
                form = "%d-%b-%Y (%H:%M:%S.%f)"
                datRef = datetime.datetime.strptime(FireSample['reference'], form)
                todayDiff = currentTime - datRef
                # Get the difference in hours
                days, seconds = todayDiff.days, todayDiff.seconds
                todayRef = days * 24 + seconds / 3600
                # Check the number of missing days
                #if days - int(FireSample['lastTimeRef']) > 1:
                #    FireSample['missing'] = FireSample['missing'] + 1
                dist = int(days - int(FireSample['lastTimeRef']) - int(FireSample['missing']))
                if dist > 1:
                    FireSample['missing'] = FireSample['missing'] + dist - 1
                # Determine the day index to append
                index = str(days - FireSample['missing'])
                # Check if there is information on the index
                if FireSample['react'].get(index) == None:
                    # Determine initial value vector
                    if currentTime.hour <= 6:
                        FireSample['initial'].append(todayRef)
                    else:
                        FireSample['initial'].append(todayRef - 1.5)
                    # Determine the sleeping moment
                    lastInd = str(int(index) - 1)
                    FireSample['final'].append(FireSample['time'][lastInd][-1])
                    # Include into react, time and error
                    FireSample['react'][index]= [float(FireLast['react'])]
                    FireSample['time'][index]= [todayRef]
                    FireSample['error'][index] = [float(FireLast['error'])]
                else:
                    # Include into react, time and error
                    FireSample['react'][index].append(float(FireLast['react']))
                    FireSample['time'][index].append(todayRef)
                    FireSample['error'][index].append(float(FireLast['error']))
                # Json dumps the lists
                FireInfo = {'time': json.dumps(FireSample['time']),
                            'react': json.dumps(FireSample['react']),
                            'error': json.dumps(FireSample['error']),
                            'initial': json.dumps(FireSample['initial']),
                            'final': json.dumps(FireSample['final']),
                            'reference': JSONsamples['reference'],
                            'lastTimeRef': index,
                            'missing': FireSample['missing'],
                            'measures': float(FireSample['measures']) + 1}
            else:
                # Initializing the first reference at 04:00:00 AM
                timeRef = currentTime.replace(hour = 4)
                timeRef = timeRef.replace(minute = 0)
                timeRef = timeRef.replace(second = 0)
                timeRef = timeRef.replace(microsecond = 0)
                # Determine the proper time reference in string
                dateInStr = timeRef.strftime("%d-%b-%Y (%H:%M:%S.%f)")
                # Determine the relative initial time
                diff = currentTime - timeRef
                # Get the difference in hours
                days, seconds = diff.days, diff.seconds
                initialRef = days * 24 + seconds / 3600
                # Create the initial time and alert dict
                __react = {"0": [float(FireLast['react'])]}
                __error = {"0": [float(FireLast['error'])]}
                __time = {"0":[initialRef]}
                FireInfo = {'time': json.dumps(__time),
                            'reference': dateInStr,
                            'react': json.dumps(__react),
                            'error': json.dumps(__error),
                            'initial': json.dumps([initialRef]),
                            'final': json.dumps([]), 'missing': 0,
                            'measures': 1, 'lastTimeRef': 0}
            fireDB.child("users").child(user_code).child("pvt").set(FireInfo)
        return FireInfo

class fireBaseEstimate(Resource):
    def get(self):
        requested_data = request.params
        # Retrieve user and algorithm information
        user_code = requested_data['user']
        algorithm = requested_data['algorithm']
        # Read the data from firebase
        JSONsamples = fireDB.child("users").child(user_code).child("kss").get().val()
        # Check if there is data on the server
        checkNull = (JSONsamples == None)
        if not(checkNull):
            count = JSONsamples['measures']
            if (algorithm == "Annealing") and count > 10 :
                # Retrive data from FireBase
                samples = {'time': json.loads(JSONsamples['time']),
                           'levAl': json.loads(JSONsamples['alert']),
                           'initial': json.loads(JSONsamples['initial']),
                           'final': json.loads(JSONsamples['final'])}
                # Estimate the model
                modelAnnealing.fit(samples)
                # Determine the parameters
                parInfo = {'omega': modelAnnealing.omega,
                           'tau': modelAnnealing.tau,
                           'phi': modelAnnealing.phi,
                           'M': modelAnnealing.M,
                           'DC': modelAnnealing.DC,
                           'y_': modelAnnealing.y_,
                           'tau_e': modelAnnealing.tau_e,
                           'algorithm': algorithm,
                           'auth': True}
                # Send data to Firebase
                fireDB.child("parameters").child(user_code).set(parInfo)
                # Send debug Info
                return "479: Sucess"
            elif (algorithm == "MOLI") and count > 10:
                return "480: Not ready."
            else:
                return "481: Not enouth data or model not valid."
        else:
            return "482: There is no K.S.S. data."

class fireBaseValidate(Resource):
    def get(self):
        requested_data = request.params
        # Retrieve user and algorithm information
        user_code = requested_data['user']
        #
        JSONpar = fireDB.child("parameters").child(user_code).get().val()
        # Determine the model parameters
        modelAnnealing.omega = float(JSONpar['omega'])
        modelAnnealing.tau = float(JSONpar['tau'])
        modelAnnealing.phi = float(JSONpar['phi'])
        modelAnnealing.M = float(JSONpar['M'])
        modelAnnealing.DC = float(JSONpar['DC'])
        modelAnnealing.y_ = float(JSONpar['y_'])
        modelAnnealing.tau_e = float(JSONpar['tau_e'])
        # Read the validation data from fireBase
        JSONsamples = fireDB.child("users").child(user_code).child("kss").get().val()
        # Generate the initial and final windows
        initial = json.loads(JSONsamples['initial'])
        final = json.loads(JSONsamples['final'])
        index = min(len(initial), len(final))
        samples = {'initial': initial[:index],
                   'final': final[:index]}
        firstAlert = json.loads(JSONsamples['alert'])["0"][0]
        # Simulate the model
        data = modelAnnealing.simulate(samples, firstAlert)
        # Transform the data to be saved
        fireInfo = {"dalert": json.dumps(data['dAl']),
                    "dtime": json.dumps(data['dtime']),
                    "nalert": json.dumps(data['nAl']),
                    "ntime": json.dumps(data['ntime'])}
        # Save data at firebase
        #fireDB.child("simulated").child(user_code).set(fireInfo)
        return fireInfo

class fireBasePredict(Resource):
    def get(self):
        requested_data = request.params
        # Retrieve the usefull info for the algorithm
        user_code = requested_data['user']
        daysAhead = requested_data['days']
        print("eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee" + daysAhead)
        alert = float(requested_data['alert'])
        # Determine the current time
        currentTime = datetime.datetime.now()
        currentTime = currentTime - datetime.timedelta(hours=7)
        # Determine the time reference from fireBase
        JSONsamples = fireDB.child("users").child(user_code).child("kss").get().val()
        referenceInStr = JSONsamples['reference']
        form = "%d-%b-%Y (%H:%M:%S.%f)"
        datRef = datetime.datetime.strptime(referenceInStr, form)
        todayDiff = currentTime - datRef
        # Get the difference in hours
        days, seconds = todayDiff.days, todayDiff.seconds
        # Determine the moment in hours
        currentHours =  days * 24 + seconds / 3600
        initial = []
        final = []
        for i in range(int(daysAhead)):
            if i == 0:
                initial.append(currentHours)
                final.append(initial[i] + 16 - currentHours % 24)
            else:
                initial.append(final[i-1] + 8)
                final.append(initial[i] + 16)
        # Read the user parameters
        JSONpar = fireDB.child("parameters").child(user_code).get().val()
        # Determine the model parameters
        modelAnnealing.omega = float(JSONpar['omega'])
        modelAnnealing.tau = float(JSONpar['tau'])
        modelAnnealing.phi = float(JSONpar['phi'])
        modelAnnealing.M = float(JSONpar['M'])
        modelAnnealing.DC = float(JSONpar['DC'])
        modelAnnealing.y_ = float(JSONpar['y_'])
        modelAnnealing.tau_e = float(JSONpar['tau_e'])
        #
        samples = {"initial": initial,
                   "final": final}
        print(samples)
        initalert = 13*(10 - float(alert))/9 + 1
        # Simulate the system
        data = modelAnnealing.simulate(samples, initalert)
        # Transform the data to be saved
        fireInfo = {"dalert": json.dumps(data['dAl']),
                    "dtime": json.dumps(data['dtime']),
                    "nalert": json.dumps(data['nAl']),
                    "ntime": json.dumps(data['ntime'])}
        return fireInfo


# Apend API resource
#api.add_resource(ToDoSimple, '/MOLI/<string:model_id>')
#api.add_resource(resAnnealing, '/Annealing/<string:model_id>')
api.add_resource(fireBaseValidate, '/validate/')
api.add_resource(fireBaseEstimate, '/estimate/')
api.add_resource(fireBaseUpdate, '/update/')
api.add_resource(fireBasePredict, '/predict/')

if __name__ == '__main__':
    app.run(port=5000, debug = True)
