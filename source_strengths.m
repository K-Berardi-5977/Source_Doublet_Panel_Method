function sigma = source_strengths(body,flow)

sigma = -(body.normal * flow.Uvec.');

end