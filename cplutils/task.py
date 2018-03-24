import time

class Info():
    def __init__(self):
        pass

class Task:
    FAILURE = 0
    SUCCESS = 1
    RUNNING = 2
    def __init__(self):
        self.compute_time = 0.0
        self.status = None
        self.info = Info()
        self.print_info = False
        self.print_progress = False
        self.print_every = 10
        self.progress = 0.0
        self.capture_exceptions = False
    
    def _run(self):
        pass

    def run(self):
        tini = time.time()
        self.status = Task.RUNNING
        self._run()
        try:
            self._run()
            self.status = Task.SUCCESS
            tfin = time.time()
            self.compute_time = tfin - tini
        except Exception as e:
            raise e
            self.status = Task.FAILURE
            if not self.capture_exceptions:
                raise e


