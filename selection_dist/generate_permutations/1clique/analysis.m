A = load("vectors.txt");
[coeff,score,latent] = princomp(A);

threshold = latent(1)*0.5;

count = [ size(find(latent > threshold))(1) size(find(latent > mean(latent)))(1) ];

count

mu=mean(latent);

j=find(latent>0);
latent = log(latent(j)/mu);
latent = latent(find(latent>0));

i=1:size(latent)(1);

weightedlatent = latent(i).*i';

2 * sum(weightedlatent) / sum(latent)

save('-ascii', 'eigenvalues.txt', 'latent');