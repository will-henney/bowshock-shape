import numpy as np

class StandingWave(object):

    def __init__(self, base_shape, amplitude, wavenumber, phase=0.0, *base_shape_args):
        """
        Constant amplitude standing wave with antinode on x-axis (theta = 0)

        Parameters
        ----------
        base_shape : callable
            Underlying R(theta) shape to be perturbed
        amplitude
            Relative amplitude of radial oscillation
        wavenumber
            Angular wavenumber of oscillation. With `wavenumber = 1.0` there will
            be a node on the y-axis (theta = pi/2). With `wavenumber = 2.0` there
            will be an anti-node on the y-axis that is out of phase with the
            node on the x-axis.
        phase : optional
            Temporal phase of the oscillation. The amplitude is multiplied
            by `cos(2 pi phase)`, so that When `phase = 0` the radial perturbation
            is positive at `theta = 0`, whereas when `phase = 0.5` the radial
            perturbation is negative at `theta = 0`.  For `phase = 0.25, 0.75`
            the perturbation is zero for all `theta`.
        *base_shape_args : optional
            Any unrecognised args are passed on to the unperturbed radius
            function: `base_shape(theta, *base_shape_args)` 


        """
        self.base_shape = base_shape
        self.amplitude = amplitude
        self.wavenumber = wavenumber
        self.phase = phase
        self.base_shape_args = base_shape_args

    def __call__(self, theta):
        """
        Radius as function of `theta` for standing wave
        """
        return ((1.0 + self.perturbation(theta)) *
                self.base_shape(theta, *self.base_shape_args))

    def perturbation(self, theta):
        """
        Fractional perturbation of `base_shape` radius at each angle `theta`
        """
        return (self.amplitude *
                np.cos(2*np.pi*self.phase) *
                np.cos(self.wavenumber*theta))
