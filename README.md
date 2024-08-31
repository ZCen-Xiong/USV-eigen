# indirect-finite-thrust-bang-bang-control
## Dynamics
The state and control variables are defined as:
$$
X = [x_1, x_2] = [R, \alpha, v, q]
$$

$$
U = [a, \omega]
$$

The state equation is given by the following:

$$
\begin{cases}
\dot{x}_1 = [ \cos q \cdot v, \sin q \cdot v] / R \\
\dot{x}_2 = [a, \omega]
\end{cases}
$$

## indirect time-optimal control
In the optimization problem setting, the state variables are set to be the following, with the state variables \( X_1 = [x_{11}, x_{12}] \), and the control variables \( X_2 = [x_{21}, x_{22}] \), which are considered for optimization.

The control input \( u \) is equal to \( X_2 \), and a feedback control law is designed to ensure the system's stability and to achieve the desired dynamic performance.

The Hamiltonian function, which is necessary for solving the optimization problem, is defined as follows:

$$
H = 1 + \lambda_1^T X_1 + \lambda_2^T X_2 = 1 + \lambda_1^T f(X_1, X_2) + \lambda_2^T g(u)
$$

The optimal control strategy is determined by the following equations:

$$
\dot{\lambda}_1 = - \frac{\partial H}{\partial X_1} 
$$

and

$$
\dot{\lambda}_2 = - \frac{\partial H}{\partial X_2}  
$$ 



The control law is designed to switch based on the value of \( \lambda_2 \) and a predefined threshold \( \varepsilon \), as follows:

$$
u = \begin{cases}
u_b, & \lambda_2 > \varepsilon \\
u_b, & \lambda_2 < \varepsilon \\
0, & \text{otherwise}
\end{cases}
$$

The system is then subjected to various simulations to verify the effectiveness of the proposed control strategy.
