function links = CartLinks(sol,setup)

 left_phase = sol.left.phase;
 right_phase = sol.right.phase;

 state_left = sol.left.state;
 state_right = sol.right.state;

 links = [state_left - state_right];
