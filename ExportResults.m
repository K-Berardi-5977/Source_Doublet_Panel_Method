
function result = ExportResults(mu, aero, mesh, numPan)

    X_Cp = mesh.Bp2(1:numPan,1);  % for plotting x locations

    result.mu   = mu;
    result.Cp   = aero.Cp;
    result.cl   = aero.cl;
    result.X_Cp = X_Cp;
end