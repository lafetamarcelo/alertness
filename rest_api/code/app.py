from flask import Flask, request
from flask_restful import Resource, Api
from matchable_observable import MOLIwhiteBox
import json

app = Flask(__name__)
api = Api(app)

class ToDoSimple(Resource):
    def get(self, model_id):
        samples = {'time': request.form['time'],
                   'levAl': request.form['levAl'],
                   'windDecisions': request.form['windDecisions']}
        return samples

    def put(self, model_id):
        samples = {'time': json.loads(request.form['time']),
                   'levAl': json.loads(request.form['levAl']),
                   'windDecisions': json.loads(request.form['windDecisions'])}
        modelMOLI = MOLIwhiteBox('trivial')
        modelMOLI.fit(samples)
        return {'omega': modelMOLI.tau}
        #toDos[toDo_id] = request.form['data']
        #return {toDo_id: toDos[toDo_id]}

api.add_resource(ToDoSimple, '/<string:model_id>')

if __name__ == '__main__':
    app.run(port=5000, debug = True)
