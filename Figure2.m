%% Code to reproduce figure 2
% Dependencies include generalized_phase.m from https://github.com/mullerlab/generalized-phase
% bandpass_filter.m from https://github.com/mullerlab/wave-matlab

load('F2_upload_data_2024update.mat')

% a - example lfp 
lfp = F2.lfp_a; 
filtLFP = bandpass_filter(reshape(lfp,1,1,[]),5,50,4,1000);

figure; hold on; plot( lfp, 'linewidth', 1, 'color', [216,216,216]/255);
plot( squeeze(filtLFP), '-k', 'linewidth',1)
xlabel( 'time (ms)'); ylabel( 'potential') 
plot([480, 480], [-2,2.5],'-k')
plot([550,550], [-2,2.5],'-k')
legend( 'raw', 'filtered')

% b - epoch
lfp = F2.lfp_b;
filtLFP = bandpass_filter(lfp,5,50,4,1000);
xgp = generalized_phase( filtLFP, 1e3,5);

chans = [28,20; 28,22; 29,19; 29,20; 29,21; 30,18; 30,19; 30,20; 30,21; 31, 18; 31,19;  31,21; 31,22; 32,18; 32,20];    
figure; hold on;
for ii = 1:length(chans), plot( squeeze( cos(angle( xgp( chans(ii,1), chans(ii,2), 480:550))) ), 'color', [216,216,216]/255); end
plot( squeeze(cos(angle( xgp(30,20, 480:550) ))),'-k', 'linewidth',1 )

% panels c/d - examples 
for ii = 1:size(F2.examples,3)
figure; imagesc( angle( F2.examples(:,:,ii)) ); map = colorcet( 'C2' );  colormap( circshift( map, [ 28, 0 ] ) ); caxis([-pi,pi]); c = colorbar; c.Label.String = 'phase (rad)';
end

% panel e - similarity matrix
figure; imagesc( F2.rho); colorbar; caxis( [-1,1])

% panel f - projection 
prj = F2.projection;
figure; scatter3( prj(:,1),prj(:,2),prj(:,3), 20, F2.colors, 'filled')

