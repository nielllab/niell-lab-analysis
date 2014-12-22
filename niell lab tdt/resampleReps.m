function [shuffParams shuffResp] = resampleReps(params,resp,reps)
trials = zeros(size(resp));
for n = 1:length(params)/reps
    trials((n-1)*reps+1 :n*reps) = (n-1)*reps + ceil(rand(reps,1)*reps);
end
% figure
% plot(trials)
shuffParams = params(trials,:);
shuffResp = resp(trials);

    