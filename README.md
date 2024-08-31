#### indirect-finite-thrust-bang-bang-control

## Requirements

This project requires the following dependencies:
- **CMake**: Minimum version required is **3.10.0**.
- **C++ Standard**: The project uses **C++17** features.
- **Eigen**: Version **3.4.0** is required for matrix and vector operations.

### Installing Eigen 3.4.0

To install Eigen 3.4.0, you can follow these steps:

#### 1. Download Eigen 3.4.0

You can download Eigen 3.4.0 from the official website or clone the repository using git:

```bash
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar -xvzf eigen-3.4.0.tar.gz
```

Alternatively, using git:

```bash
git clone https://gitlab.com/libeigen/eigen.git
cd eigen
git checkout 3.4.0
```

#### 2. Install Eigen 3.4.0

Eigen is a header-only library, so you don't need to compile it. Just include the Eigen headers in your project. To make it available globally:

```bash
cd eigen-3.4.0
mkdir build
cd build
cmake ..
sudo make install
```

This will install Eigen headers into `/usr/local/include/eigen3/`.

### CMake Configuration

Make sure your `CMakeLists.txt` is set to require the correct CMake version, C++ standard, and links the Eigen library properly:

```cmake
cmake_minimum_required(VERSION 3.10)
project(indirect-finite-thrust-bang-bang-control)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 3.4.0 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(your_executable_name main.cpp)
target_link_libraries(your_executable_name Eigen3::Eigen)
```
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
