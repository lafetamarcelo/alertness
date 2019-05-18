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
        # Check the results for Classification
        if (examsCond[0]):
            # Retrieve data samples from FireBase
            JSONsamples = fireDB.child("users").child(user_code).child("class").get().val()
            if JSONsamples != None:
                # Format the data into a dictionary
                FireSample = {'time': json.loads(JSONsamples['time']),
                              'reference': JSONsamples['reference'],
                              'classes': json.loads(JSONsamples['classes'])}
                # Determine the hours expent since reference
                form = "%d-%b-%Y (%H:%M:%S.%f)"
                datRef = datetime.datetime.strptime(FireSample['reference'], form)
                todayDiff = datetime.datetime.now() - datRef
                # Get the difference in hours
                days, seconds = todayDiff.days, todayDiff.seconds
                todayRef = days * 24 + seconds / 3600
                # Include the class into classes
                FireSample['classes'][str(days)].append(FireLast['classe'])
                # Append the info into the time serie
                FireSample['time'][str(days)].append(todayRef)
                # Json dumps the lists
                FireInfo = { 'time': json.dumps(FireSample['time']),
                            'reference': JSONsamples['reference'],
                            'classes': json.dumps(FireSample['classes'])}
            else:
                # Format the actual time in a particular reference
                dateInStr = datetime.datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")
                # Create the initial time and alert dict
                __class = {"0": [FireLast["classe"]]}
                __time = {"0": [0]}
                FireInfo = {'time': json.dumps(__time),
                              'reference': dateInStr,
                              'classes': json.dumps(__class)}
            # Set the data back in firebase
            fireDB.child("users").child(user_code).child("class").set(FireInfo)

        if (examsCond[1]):
            # Retrieve data samples from fireBase
            JSONsamples = fireDB.child("users").child(user_code).child("kss").get().val()
            # Transfotm sleepiness(score) into alertness
            scoreToAlert = 13*(10 - double(FireLast['score']))/9 + 1
            # Check if there is info on KSS for [user_code]
            if JSONsamples != None:
                # Format data as dictionary
                FireSample = {'time': json.loads(JSONsamples['time']),
                              'alert': json.loads(JSONsamples['alert']),
                              'reference': JSONsamples['reference']}
                # Determine the hours expent since reference
                form = "%d-%b-%Y (%H:%M:%S.%f)"
                datRef = datetime.datetime.strptime(FireSample['reference'], form)
                todayDiff = datetime.datetime.now() - datRef
                # Get the difference in hours
                days, seconds = todayDiff.days, todayDiff.seconds
                todayRef = days * 24 + seconds / 3600
                # Determine the day index
                ind = days
                # Include into alert
                FireSample['alert'][str(ind)].append(scoreToAlert)
                # Append the info into the time serie
                FireSample['time'][str(ind)].append(todayRef)
                # Json dumps the lists
                FireInfo = {'time': json.dumps(FireSample['time']),
                            'reference': JSONsamples['reference'],
                            'alert': json.dumps(FireSample['alert'])}
            else:
                # Format the actual time in a particular reference
                dateInStr = datetime.datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")
                # Create the initial time and alert dict
                __alert = {"0": [scoreToAlert]}
                __time = {"0":[0]}
                FireInfo = {'time': json.dumps(__time),
                            'reference': dateInStr,
                            'alert': json.dumps(__alert)}
            fireDB.child("users").child(user_code).child(kss).set(FireInfo)
            
        if (examsCond[2]):
            # Retrieve data samples from fireBase
            JSONsamples = fireDB.child("users").child(user_code).child("pvt").get().val()
            # Check if there is info on PVT for [user_code]
            if JSONsamples != None:
                # Format data as dictionary
                FireSample = {'time': json.loads(JSONsamples['time']),
                              'react': json.loads(JSONsamples['react']),
                              'reference': JSONsamples['reference'],
                              'error': json.loads(JSONsamples['error'])}
                # Determine the hours expent since reference
                form = "%d-%b-%Y (%H:%M:%S.%f)"
                datRef = datetime.datetime.strptime(FireSample['reference'], form)
                todayDiff = datetime.datetime.now() - datRef
                # Get the difference in hours
                days, seconds = todayDiff.days, todayDiff.seconds
                todayRef = days * 24 + seconds / 3600
                # Determine the day index
                ind = days
                # Include into react, time and error
                FireSample['react'][str(ind)].append(FireLast['react'])
                FireSample['time'][str(ind)].append(todayRef)
                FireSample['error'][str(ind)].append(FireLast['error'])
                # Json dumps the lists
                FireInfo = {'time': json.dumps(FireSample['time']),
                            'reference': JSONsamples['reference'],
                            'react': json.dumps(FireSample['react']),
                            'error': json.dumps(FireSample['error'])}
            else:
                # Format the actual time in a particular reference
                dateInStr = datetime.datetime.now().strftime("%d-%b-%Y (%H:%M:%S.%f)")
                # Create the initial time and alert dict
                __react = {"0": [FireLast['react']]}
                __error = {"0": [FireLast['error']]}
                __time = {"0":[0]}
                FireInfo = {'time': json.dumps(__time),
                            'reference': dateInStr,
                            'react': json.dumps(__react),
                            'error': json.dumps(__error)}
            fireDB.child("users").child(user_code).child("pvt").set(FireInfo)
        return FireInfo

# Apend API resource
#api.add_resource(ToDoSimple, '/MOLI/<string:model_id>')
#api.add_resource(resAnnealing, '/Annealing/<string:model_id>')
#api.add_resource(resEstimate, '/estimate/')
api.add_resource(fireBaseUpdate, '/update/')

if __name__ == '__main__':
    app.run(port=5000, debug = True)
