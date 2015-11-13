%%%analyzeGeometrics
%%%use to get spiking relative to geometric stimulus movie
%%% P.R.L. Parker, Cris Niell Lab, 11/12/2015

%use getMovieData to pull spike times, mouse velocity/time, frame info
[spikes mouseT mouseV framenums frametimes] = getMovieData;

%load movie legend
load GeomStim