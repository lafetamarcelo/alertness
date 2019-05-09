from flask import Flask, request, jsonify
from flask_restful import Resource, Api
from flask_request_params import bind_request_params
from matchable_observable import MOLIwhiteBox, forceEstimate
import json
import pyrebase

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

# Apend API resource
api.add_resource(ToDoSimple, '/MOLI/<string:model_id>')
api.add_resource(resAnnealing, '/Annealing/<string:model_id>')
api.add_resource(resEstimate, '/estimate/')

if __name__ == '__main__':
    app.run(port=5000, debug = True)
