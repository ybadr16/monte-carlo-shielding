from manim import *
import json
import random

class ParticleSimulation(Scene):
    def construct(self):
        # Load JSON data
        with open("all_trajectories.json", "r") as f:
            data = json.load(f)

        # Define shield bounds
        shield_x_bounds = (-10, 10)
        shield_y_bounds = (-20, 20)
        shield_z_bounds = (-10, 10)

        def within_bounds(position):
            x, y, z = position
            return (
                shield_x_bounds[0] < x <= shield_x_bounds[1]
                and shield_y_bounds[0] < y <= shield_y_bounds[1]
                #and shield_z_bounds[0] < z <= shield_z_bounds[1]
            )
        def is_absorbed(history):
            # Absorbed if the particle stays within bounds and ends inside
            return all(within_bounds(pos) for pos in history) and within_bounds(history[-1])

        def is_scattered(history):
            # Scattered if the particle leaves the bounds at any point
            return any(not within_bounds(pos) for pos in history)


        # Filter particle histories
        histories = list(data.values())
        selected_histories = []

        for history in histories:
            if is_scattered(history) or is_absorbed(history):
                selected_histories.append(history)
            if len(selected_histories) >= 100:
                break

        scale_factor = 1

        # Create the shield
        shield_width = shield_x_bounds[1] - shield_x_bounds[0]
        shield_height = shield_y_bounds[1] - shield_y_bounds[0]
        rectangle = Rectangle(
            height=shield_height / scale_factor,
            width=shield_width / scale_factor,
            color=GRAY
        )
        rectangle.set_fill(BLUE, opacity=0.2)
        self.play(Create(rectangle))

        # Animate particle motion
        for history in selected_histories:
            initial_position = history[0]
            scaled_x = initial_position[0] / scale_factor
            scaled_y = initial_position[1] / scale_factor
            scaled_z = initial_position[2] / scale_factor

            neutron = Dot(point=np.array([scaled_x, scaled_y, 0]), color=YELLOW)
            self.add(neutron)

            for i in range(len(history) - 1):
                end = history[i + 1]
                scaled_x = end[0] / scale_factor
                scaled_y = end[1] / scale_factor
                scaled_z = end[2] / scale_factor

                self.play(neutron.animate.move_to(np.array([scaled_x, scaled_y, 0])), run_time=0.2)

            self.play(FadeOut(neutron), run_time=0.2)

        self.wait()


class ParticleSimulation3D(ThreeDScene):
    def construct(self):
        # Load JSON data
        with open("all_trajectories.json", "r") as f:
            data = json.load(f)

        # Subsample particle histories
        histories = list(data.values())
        selected_histories = random.sample(histories, 50)  # Select 50 particles for performance

        # Define shield bounds
        shield_x_bounds = (-10, 10)
        shield_y_bounds = (-10, 10)
        shield_z_bounds = (-10, 10)

        # Create the shield as a 3D rectangle (thin slab)
        shield = Cube(
            side_length=shield_y_bounds[1] - shield_y_bounds[0],
            fill_color=RED,
            fill_opacity=0.2,
            stroke_color=RED,
        )
        shield.set_width((shield_x_bounds[1] - shield_x_bounds[0]) / 5)
        shield.set_depth((shield_z_bounds[1] - shield_z_bounds[0]) / 5)
        shield.move_to(np.array([0, 0, 0]))  # Centered at x = 0
        self.play(Create(shield))

        # Set up the 3D camera
        self.set_camera_orientation(phi=0 * DEGREES, theta= 0* DEGREES)

        # Animate particle motion
        for history in selected_histories:
            positions = history  # Each particle's trajectory

            # Create a dot for the particle at its initial position
            initial_position = positions[0]
            x, y, z = initial_position

            # Scale positions for visualization
            scaled_x = x / 5
            scaled_y = y / 5
            scaled_z = z / 5

            neutron = Sphere(radius=0.1, color=YELLOW).move_to(
                np.array([scaled_x, scaled_y, scaled_z])
            )
            self.add(neutron)

            # Move the particle along its trajectory
            for i in range(len(positions) - 1):
                start = positions[i]
                end = positions[i + 1]

                # Calculate displacement
                dx = (end[0] - start[0]) / 5  # Scale for visualization
                dy = (end[1] - start[1]) / 5
                dz = (end[2] - start[2]) / 5

                # Animate motion in 3D
                self.play(neutron.animate.shift(np.array([dx, dy, dz])), run_time=0.1)

            # Remove the neutron after it completes its trajectory
            self.play(FadeOut(neutron), run_time=0.2)

        self.wait()
