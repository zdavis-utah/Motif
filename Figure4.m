% Code to Reproduce Figure 4
% August 2024

load( 'F4_upload_data.mat')

% panel c 

% panel c top 
trace = F4.trace; t = F4.t; 
figure; hold on; 
for ii = 1: size(trace,1)
    plot( trace( ii, t-50:t+50) , 'linewidth', 1, 'color', [200,200,200]*(1/255))
end
plot( trace( 25, t-50:t+50) , '-k','linewidth', 1)
xlim( [1, 100]);  ylabel( 'potential (mV)') ;  xlabel( 'time (ms)')

% panel c bottom 
map = colorcet( 'C2' ); 
figure; imagesc( F4.c_phase_plot); colormap( circshift( map, [ 28, 0 ] ) )
set(gca, 'XTick', [], 'YTick', [])
xlabel( 'Channel'); ylabel( 'Channel') 
c = colorbar; c.Label.String = 'Phase (rad)';

% panel d

% repetition 1
figure; imagesc( F4.d_rep1);  colormap( circshift( map, [ 28, 0 ] ) )%colormap(C); 
set(gca, 'XTick', [], 'YTick', [])
xlabel( 'Channel'); ylabel( 'Channel') 
c = colorbar; c.Label.String = 'Phase (rad)';

% repetition 2
figure; imagesc( F4.d_rep2);  colormap( circshift( map, [ 28, 0 ] ) )%colormap(C); 
set(gca, 'XTick', [], 'YTick', [])
xlabel( 'Channel'); ylabel( 'Channel') 
c = colorbar; c.Label.String = 'Phase (rad)';

% panel e

% non repetition 1
figure; imagesc( F4.e_nonrep1);  colormap( circshift( map, [ 28, 0 ] ) )%colormap(C); 
set(gca, 'XTick', [], 'YTick', [])
xlabel( 'Channel'); ylabel( 'Channel') 
c = colorbar; c.Label.String = 'Phase (rad)';

% non repetition 2 
figure; imagesc( F4.e_nonrep2);  colormap( circshift( map, [ 28, 0 ] ) )%colormap(C); 
set(gca, 'XTick', [], 'YTick', [])
xlabel( 'Channel'); ylabel( 'Channel') 
c = colorbar; c.Label.String = 'Phase (rad)';

% panel f

figure; imagesc( F4.similarity_matrix_example);
c = colorbar; c.Label.String = 'Similarity'; caxis([-1,1])

% panel g

figure; imagesc( F4.similarity_subset); c = colorbar; c.Label.String = 'Similarity';
set(gca, 'XTick', [], 'YTick', []); caxis( [-1,1])

% panel h

prj = F4.prj; 

% PROJECTION PLOT to demonstrate clustering
figure; scatter3( prj(:,1), prj(:,2), prj(:,3), 10, [0.5,0.5,0.5], 'filled', 'markerfacealpha',0.1) ; %plot all trials
hold on;
% add trials from specific motifs
scatter3( prj(F4.motif_example1,1), prj(F4.motif_example1,2), prj(F4.motif_example1,3), 20, [0 0.4470 0.7410], 'filled') ; %plot all trials
scatter3( prj(F4.motif_example2,1), prj(F4.motif_example2,2), prj(F4.motif_example2,3), 20, [0.8500 0.3250 0.0980], 'filled') ; %plot all trials
scatter3( prj(F4.motif_example3,1), prj(F4.motif_example3,2), prj(F4.motif_example3,3), 20, [0.9290 0.6940 0.1250], 'filled') ; %plot all trials
scatter3( prj(F4.motif_example4,1), prj(F4.motif_example4,2), prj(F4.motif_example4,3), 20, [0.4660 0.6740 0.1880], 'filled') ; %plot all trials
set(gca, 'Fontsize', 14)
xlabel('Dimension 1'); ylabel('Dimension 2'); zlabel( 'Dimension 3') 



