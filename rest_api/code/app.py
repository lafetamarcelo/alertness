from flask import Flask, request
from flask_restful import Resource, Api
from matchable_observable import MOLIwhiteBox

app = Flask(__name__)
api = Api(app)

toDos = {}

class ToDoSimple(Resource):
    def get(self, toDo_id):
        send = int(toDos[toDo_id]) * 2
        return {toDo_id: send}

    def put(self, toDo_id):
        toDos[toDo_id] = request.form['data']
        return {toDo_id: toDos[toDo_id]}

api.add_resource(ToDoSimple, '/toDos/<string:toDo_id>')

if __name__ == '__main__':
    app.run(port=5000, debug = True)
