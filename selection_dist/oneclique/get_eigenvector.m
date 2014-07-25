arglist = argv();

A = load(arglist{1});

[V,E] = eig(A);

E=diag(E);

E=E.*conj(E);

[E,i] = sort(E, 'descend');

V=V(:,i);

firstvec = real(V(:,1));
secondvec = real(V(:,2));

firstvecI = imag(V(:,1));
secondvecI = imag(V(:,2));

save("-ascii", "first.txt", "firstvec");
save("-ascii", "second.txt", "secondvec");
save("-ascii", "firstI.txt", "firstvecI");
save("-ascii", "secondI.txt", "secondvecI");
