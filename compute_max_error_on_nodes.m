function result=compute_max_error_on_nodes(u,uh)

err = u-uh;
abs_err = abs(err);

result = max(abs_err);
end