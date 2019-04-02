from flask import Flask, request, jsonify
from flask_restful import Resource, Api
from flask_request_params import bind_request_params
from matchable_observable import MOLIwhiteBox
import json

app = Flask(__name__)
api = Api(app)
app.before_request(bind_request_params)

modelMOLI = MOLIwhiteBox('trivial')

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
        #modelMOLI = MOLIwhiteBox('trivial')
        modelMOLI.fit(samples)
        #self.model = modelMOLI
        return {'omega': modelMOLI.w}

api.add_resource(ToDoSimple, '/<string:model_id>')

if __name__ == '__main__':
    app.run(port=5000, debug = True)
