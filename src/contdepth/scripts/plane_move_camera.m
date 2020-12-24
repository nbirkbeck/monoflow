% The shell scripts used to condense the data.
%(for i in 2 3 4 6 10; do dir=/tmp/plane_$i; echo ${dir}; ./compare --flow=${dir}/results/sf-flow-%d.rfi ${dir}/seq/{0,1,2}/image-0000.png --depth=${dir}/results/sf-depth-%d.rfi; ./compare ${dir}/seq/{0,1,2}/image-0000.png --depth=${dir}/results/vari-depth-%d.rfi; done) > /tmp/plane_move_results.txt
%grep "disp:" /tmp/plane_move_results.txt | perl -p -i -e ['s/disp:\[//g'] | gawk -F , '{if ((NR % 9) % 2 == 0 && (NR % 9) <= 6) print $1; if ((NR) % 9 == 0) print "---";}' > /tmp/plane_flow.txt
%grep "depth:" /tmp/plane_move_results.txt | perl -p -i -e 's/depth:\[//g' | gawk -F , '{print $1; if ((NR) % 6 == 0) print "---";}' > /tmp/plane_depth.txt

% Camera distances:  1.9, 1.3,  
% 4 units from object, object square 2x2, disp 0.2 normal to plane,
% 0.1 translational

d = [0.0204012, 0.0217432, 0.0322414;
     0.0111906, 0.0127422, 0.0181849;
     0.00720615, 0.00884378, 0.0107525;
     0.00472139, 0.00572837, 0.00610781;
     0.00752661, 0.00685678, 0.00732034];

v = [0.25244, 0.232245, 0.259176;
     0.136945, 0.12416, 0.138708;
     0.0948975, 0.0843385, 0.0952558;
     0.0603554, 0.0515657,0.0598086;
     0.0390689, 0.0300954, 0.0482355];


f = [0.0053678, 0.00429832, 0.0057561;
     0.00278693, 0.00254966, 0.00324918;
  0.00200743, 0.00202942, 0.00220582;
  0.00156157,0.00173133,0.00167063;
  0.00292935, 0.00296983, 0.0028755];


davg = mean(d,2);
favg = mean(f,2);
vavg = mean(v,2);

clf;
subplot(2, 2, 1);

x=(([2,3,4,6,10]-1)/9);
plot(x, vavg, 'r-.');hold on;
plot(x, davg, 'b-');
plot(x, favg * 10, 'k:o');

legend('Init (Vari)', 'Disp',  'Flow (* 10)');
xlabel('Baseline Proportion')
ylabel('Distance')