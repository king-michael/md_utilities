"""
Different fall back solutions for dependencies
"""

# Fall back for ProgressReporter_
# try:
#     from progress_reporter import ProgressReporter_
# except:
#     class ProgressReporter_:
#       ...

class ProgressReporter_:
    def __init__(self):
        """
        Examples
        --------
        >>> n_iterations = 100
        >>> pgr = ProgressReporter_()
        >>> pgr.register(n_iterations, description="test")
        >>> for i in range(n_iterations):
        >>>     pgr.update(1)
        >>> pgr.finish()
        """
        pass
    def register(self,amount_of_work, description=None, stage=0):
        self.count = 0
        self.amount_of_work = amount_of_work
        lsteps = len(str(amount_of_work))+1
        self.description = description
        if description is None:
            self.print_format = '\r{{step:{lsteps}d}}/{amount_of_work} ({{percent:5.2f}}%)'.format(**locals())
        else:
            self.print_format = '\r{description}:{lsteps}d}}/{amount_of_work} ({{percent:5.2f}}%)'.format(**locals())

    def update(self, increment, stage=0):
        self.count += increment
        print(self.print_format.format(step=self.count,
                                       percent=100.0*self.count/self.amount_of_work),
              end='')
    def finish(self):
        print('\nfinished')