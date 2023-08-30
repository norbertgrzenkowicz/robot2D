# robot2D

This repository contains a Python project developed to visualize the animation of a 2D manipulator. The manipulator's motion is simulated using techniques such as Inverse Kinematics, Linear Kinematics, and Lagrangian Dynamics. This project was created to fulfill the requirements of a university course and to provide a hands-on understanding of these fundamental concepts in robotics and mechanics.

## Project Overview

The primary objective of this project is to provide a visual representation of the motion and behavior of a 2D manipulator. It achieves this through the implementation of various techniques:

- **Inverse Kinematics**: The program calculates the joint angles required to place the manipulator's end effector at a specified position in space. This is particularly useful for tasks involving precise positioning.

- **Linear Kinematics**: Linear motion of the manipulator's end effector is simulated based on given joint velocities. This helps in understanding how joint velocities influence the overall motion of the manipulator.

- **Lagrangian Dynamics**: The manipulator's motion is simulated by solving the Lagrange equations, which involve kinetic and potential energy terms. This provides insights into the dynamic behavior of the manipulator and its response to external forces.

## Usage

To run the animation and interact with the manipulator simulation, follow these steps:

1. Ensure you have Python installed on your system.
2. Clone this repository to your local machine.
3. Navigate to the repository's directory using the terminal.
4. Run the main Python script to start the animation:

   ```bash
   python main.py
   ```

5. The animation window will open, displaying the manipulator's motion based on the implemented techniques. Interact with the animation as needed to observe different scenarios.

## Repository Structure

The repository is organized as follows:

- `main.py`: The main Python script that orchestrates the animation and simulation.
- `inverse_kinematics.py`: Implementation of inverse kinematics calculations.
- `linear_kinematics.py`: Implementation of linear kinematics simulation.
- `lagrangian_dynamics.py`: Implementation of Lagrangian dynamics simulation.
- `utils.py`: Utility functions used across the project.
- `assets/`: Directory containing any additional assets, such as images or data files.

## Acknowledgments

This project was developed as a part of a university course to provide practical insights into the concepts of inverse kinematics, linear kinematics, and Lagrangian dynamics. The animation and simulations aim to enhance understanding and intuition regarding the behavior of 2D manipulators in different scenarios.

## Note

While this project focuses on demonstrating specific concepts, please keep in mind that real-world manipulator behaviors can be more complex. This project provides a simplified representation for educational purposes.

For any questions, issues, or improvements, feel free to create an issue or pull request in the repository. Happy learning!
