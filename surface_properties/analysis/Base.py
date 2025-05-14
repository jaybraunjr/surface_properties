from tqdm import tqdm

class AnalysisBase:
    def __init__(self, universe, start=0, stop=None, step=1, verbose=True):
        self.u = universe
        self.start = start
        self.stop = stop or len(universe.trajectory)
        self.step = step
        self.verbose = verbose 
        self.results = {}
        self._frame_index = 0
        
    def _prepare(self):
        """Set up data structures. Override in subclass."""
        pass

    def _analyze_frame(self, ts):
        """Do per-frame analysis. Override in subclass."""
        raise NotImplementedError

    def _finalize(self):
        """Any cleanup or final stats. Optional."""
        pass

    def run(self):

        """Run the analysis."""
        self._prepare()

        traj = self.u.trajectory[self.start:self.stop:self.step]
        iterator = tqdm(traj, desc=self.__class__.__name__) if self.verbose else traj

        for ts in iterator:
            self._analyze_frame(ts)
            self._frame_index += 1

        self._finalize()
        return self
