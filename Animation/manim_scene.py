from manim import *
import json

class ParticleSimulation(Scene):
    def construct(self):
        # Load JSON data
        with open("all_trajectories.json", "r") as f:
            data = json.load(f)

        # Define shield bounds
        shield_x_bounds = (-10, 10)
        shield_y_bounds = (-20, 20)
        z_max = 20

        def within_bounds(position):
            x, y, z = position
            return (
                shield_x_bounds[0] < x <= shield_x_bounds[1]
                and shield_y_bounds[0] < y <= shield_y_bounds[1]
                and z <= z_max
            )

        def is_valid_particle(history):
            # Exclude particles that exit due to z > 20
            if any(pos[2] > z_max for pos in history):
                return False

            # Adjust the last coordinate for particles with more than two entries
            #if len(history) > 2:
                # Check if the last coordinate exits the shield's x or y bounds
            #    if not within_bounds(history[-1]):
                    # Use the second-to-last coordinate as the final position
            #        history[-1] = history[-2]

            return True

        # Filter particle histories
        selected_histories = {}
        for key, history in data.items():
            if is_valid_particle(history):
                selected_histories[key] = history
            if len(selected_histories) >= 100:
                break

        scale_factor = 7

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

        # Add text for displaying the current particle key
        particle_key_text = Text("", font_size=24).to_corner(UL)
        self.add(particle_key_text)

        # Animate particle motion
        for key, history in selected_histories.items():
            # Update the displayed key
            particle_key_text.become(Text(f"Particle Key: {key}", font_size=24).to_corner(UL))

            initial_position = history[0]
            scaled_x = initial_position[0] / scale_factor
            scaled_y = initial_position[1] / scale_factor

            neutron = Dot(point=np.array([scaled_x, scaled_y, 0]), color=YELLOW)
            self.add(neutron)

            for i in range(len(history) - 1):
                end = history[i + 1]
                scaled_x = end[0] / scale_factor
                scaled_y = end[1] / scale_factor

                self.play(neutron.animate.move_to(np.array([scaled_x, scaled_y, 0])), run_time=0.2)

            self.play(FadeOut(neutron), run_time=0.2)

        self.wait()
