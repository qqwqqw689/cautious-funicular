#Pendulum simulation
import numpy as np
import matplotlib.pyplot as plt


class ODESolver(object):
	"""Second_order OED Solver.
		Parameters
    	------------
    	omega_0 : float
    	        initial angular velocity
    	theta_0 : float
    	        initial angular displacement
    	eta : float
    	    time step size
    	n_iter : int
    	       number of steps
	
    	Attributes
    	-----------
    	time_ : 1d-array
    	    Stores time values for each time step.
    	omega_ : 1d-array
    	    Stores angular velocity values for each time step.
    	theta_ : 1d-arra
    	   Stores angular displacement values for each time step.
	
    	Methods
    	-----------
    	euler(alpha): Implements the Euler algorithm for the acceleration function alpha.
    
    	midpoint(alpha): Implements the Midpoint algorithm for the acceleration function alpha.
    
    	verlet(alpha): Implements the Verlet algorithm for the acceleration function alpha.
    """
	def __init__(self,omega_0=0,theta_0=10,eta=0.01,n_iter=10):
		self.omega_0=omega_0
		self.theta_0=theta_0
		self.eta=eta
		self.n_iter=n_iter

	def euler(self,alpha):
		"""
		Implements Euler Method.
    
    	Parameters
    	----------
    	alpha : acceleration function

    	Returns
    	-------
    	self : object
    	"""
		self.time_=np.zeros(self.n_iter)
		self.omega_=np.zeros(self.n_iter)
		self.theta_=np.zeros(self.n_iter)
		self.omega_[0]=self.omega_0
		self.theta_[0]=self.theta_0*np.pi/180.0
		for i in range(self.n_iter-1):
			self.time_[i+1]=self.time_[i]+self.eta
			self.omega_[i+1]=self.omega_[i]+self.eta*alpha(self.theta_[i])
			self.theta_[i+1]=self.theta_[i]+self.eta*self.omega_[i]
		
		return self

	def midpoint(self,alpha):
		"""
		Implements Midpoint Method.
    
    	Parameters
    	----------
    	alpha : acceleration function

    	Returns
    	-------
    	self : object
    	"""
		self.time_=np.zeros(self.n_iter)
		self.omega_=np.zeros(self.n_iter)
		self.theta_=np.zeros(self.n_iter)
		self.omega_[0]=self.omega_0
		self.theta_[0]=self.theta_0*np.pi/180.0
		for i in range(self.n_iter-1):
			self.time_[i+1]=self.time_[i]+self.eta
			self.omega_[i+1]=self.omega_[i]+self.eta*alpha(self.theta_[i])
			self.theta_[i+1]=self.theta_[i]+0.5*self.eta*(self.omega_[i]+self.omega_[i+1])
		return self


#Define angular Acceleration
def alpha(x):
	return -np.sin(x)

time=ODESolver(omega_0=0,theta_0=10,eta=0.000001,n_iter=30000000).euler(alpha).time_
theta=ODESolver(omega_0=0,theta_0=10,eta=0.000001,n_iter=30000000).euler(alpha).theta_

plt.plot(time,theta*180/np.pi,lw=3,color='red')

#   lw : Set the line width in points
plt.xlabel('time(s)')
plt.ylabel('angle(deg)')
plt.title('Euler Method')
plt.show()


time=ODESolver(omega_0 = 0, theta_0 = 10, eta=0.1, n_iter=300).midpoint(alpha).time_
theta=ODESolver(omega_0 = 0, theta_0 = 10, eta=0.1, n_iter=300).midpoint(alpha).theta_
plt.plot(time,theta*180/np.pi,lw=3,color='green')
plt.xlabel('time(s)',size=13)
plt.ylabel('angle (deg)',size=13)
plt.title('Midpoint Method',size=13)
plt.show()
