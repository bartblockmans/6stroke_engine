function [] = display_progress(phase)
% phase = 1 : load parameters
% phase = 2 : precompute engine parameters & geometry
% phase = 3 : integrate engine EOMs & find periodic steady state
% phase = 4 : include effect of DMF & compute gearbox input torque
% phase = 5 : visualize results
% phase = 6 : animate results

switch phase

    case 1

disp('===================================================================')
disp('           PORSCHE 6-STROKE COMBUSTION ENGINE SIMULATION           ')
disp('===================================================================')
disp(' ');
disp('1. Load simulation parameters')
disp('-------------------------------------------------------------------')
disp(' ');
disp('   Loading simulation parameters...')

    case 2

disp('   Simulation parameters loaded.')
disp(' ');
disp('2. Precompute engine parameters & geometry')
disp('-------------------------------------------------------------------')
disp(' ');
disp('   Precomputing engine parameters...')

    case 3

disp('   Engine parameters computed. Ready for simulation.')
disp(' ');
disp('3. Integrate engine EOMs and find periodic steady-state')
disp('-------------------------------------------------------------------')
disp(' ');
disp('   Fixed-point iteration on integrated engine EOMs...')
disp(' ');

    case 4

disp(' ');
disp('   Periodic steady-state of the engine computed.')
disp(' ');
disp('   Simulating complete 6-cylinder engine...')

    case 5

disp('   Engine response computed.')
disp(' ');
disp('   Plotting key simulation results...')

    case 6

disp('   Key simulation plots generated.')
disp(' ');
disp('4. Animate the simulation results')
disp('-------------------------------------------------------------------')
disp(' ')
disp('   Animation of the simulation results now shown.')
disp(' ')

end
