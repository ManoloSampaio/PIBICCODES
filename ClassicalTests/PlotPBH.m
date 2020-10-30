function PlotPBH(A,B)

num_of_states = length(A);
xticks_label = [1:1:num_of_states];
[eigen,ranks]=PBHc(A,B);
stem(ranks,'BaseValue',num_of_states,'LineStyle','none')
xlabel("Eigenvalue of the test");
ylabel("Rank");
xticks(xticks_label);
end