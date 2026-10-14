"""Compatibility facade for the v2 Trajectory/RobustnessComputer API."""
from sinabro.engine import Trajectory as _Trajectory, RobustnessComputer
from sinabro.core import MutationRecord

class Trajectory(_Trajectory):

    def __init__(self, traj_id, record, note='.', **kwargs):
        if not isinstance(record, MutationRecord):
            raise TypeError('record must be MutationRecord')
        super().__init__(traj_id, record.sequence, note=record.note if note == '.' else note, **kwargs)

    def autofill(self, method, eval_method, **kwargs):
        """Generate substitutions until a criterion, fixed length, no event, or safety cap."""
        if 'maxlen' in kwargs:
            kwargs['max_length'] = kwargs.pop('maxlen')
        condition = 'evaluation' if eval_method == 'blosum' else eval_method
        return super().autofill(condition, method=method, **kwargs)
