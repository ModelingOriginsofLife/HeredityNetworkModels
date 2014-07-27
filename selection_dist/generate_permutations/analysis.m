A = load("vectors.txt");
[coeff,score,latent] = princomp(A);

threshold = latent(1)*0.5;

j=find(latent>0);
latent = latent(j);

count = [ size(find(latent > threshold))(1) size(find(latent > mean(latent)))(1) size(find(log(latent) > mean(log(latent))))(1) ];

count

mu=mean(latent);

latent = log(latent/mu);
latent = latent(find(latent>0));

i=1:size(latent)(1);

weightedlatent = latent(i).*i';

2 * sum(weightedlatent) / sum(latent)

save('-ascii', 'eigenvalues.txt', 'latent');