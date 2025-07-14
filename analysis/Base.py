"""Core analysis base classes.

This module defines :class:`AnalysisBase`, a lightweight framework used by
all analysis routines in this repository.  It provides a standard ``run``
method which iterates over a trajectory and delegates work to three hooks:

``_prepare``
    Called once before iteration.  Subclasses typically allocate any
    arrays or perform selections here.
``_analyze_frame``
    Invoked for each timestep.  The current :class:`~MDAnalysis.core.Timestep`
    is supplied allowing subclasses to accumulate results frame by frame.
``_finalize``
    Executed after all frames have been processed to compute any aggregated
    quantities.
"""

from tqdm import tqdm


class AnalysisBase:
    """Base class for frame-by-frame analyses.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory to analyse.
    start, stop, step : int, optional
        Frame range and stride used when iterating over the trajectory.
    verbose : bool, optional
        If ``True`` a progress bar is displayed while running.

    Notes
    -----
    Subclasses are expected to override :meth:`_prepare`, :meth:`_analyze_frame`
    and :meth:`_finalize`.  Each of these methods may access or populate the
    ``results`` dictionary attribute as required.
    """
    def __init__(self, universe, start=0, stop=None, step=1, verbose=True):
        self.u = universe
        self.start = start
        self.stop = stop or len(universe.trajectory)
        self.step = step
        self.verbose = verbose 
        self.results = {}
        self._frame_index = 0
        
    def _prepare(self):
        """Set up analysis-specific state prior to iteration.

        Subclasses typically override this to create selections or allocate
        arrays that will store per-frame results.
        """
        pass

    def _analyze_frame(self, ts):
        """Perform calculations for a single frame.

        Parameters
        ----------
        ts : MDAnalysis.core.groups.Timestep
            Timestep object representing the current simulation frame.
        """
        raise NotImplementedError

    def _finalize(self):
        """Finalize the analysis after all frames have been processed.

        This hook can be used to perform averaging or to compute summary
        statistics from values accumulated during :meth:`_analyze_frame`.
        """
        pass

    def run(self):
        """Execute the analysis workflow.

        Returns
        -------
        AnalysisBase
            The instance itself to allow for fluent usage.
        """
        self._prepare()

        traj = self.u.trajectory[self.start:self.stop:self.step]
        iterator = tqdm(traj, desc=self.__class__.__name__) if self.verbose else traj

        for ts in iterator:
            self._analyze_frame(ts)
            self._frame_index += 1

        self._finalize()
        return self
    
