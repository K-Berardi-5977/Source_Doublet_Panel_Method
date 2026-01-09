function [sigma, Nuemann_Validation] = FixSourceStrength(U_normal, S)
%Simple function which fixes the values of the source strengths for each
%panel by setting them to the local normal component of the free-stream
%velocity - the function features a validation check that the Nuemann
%Boundary condition is satisfied

sigma = U_normal; % Set source strength equal to normal velocity at each panel

%Validate that the sum of all source strengths is approximately zero such
%that no net fluid is added or removed from the flow

if sum(sigma(:).*S(:)) <= 1E-4
    disp(sum(sigma.*S))
    Nuemann_Validation = true;
    disp("The Nuemann Boundary Condition is Satisfied")
else
    disp(sigma(:).*S(:))
    Nuemann_Validation = false;
    disp("The Nuemann Boundary Condition is Violated")
end

end