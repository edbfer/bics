i = 0;
while((1)){
  splot "dados/t".(i % 201).".txt" using 1:2:3 with dots;
  pause 0.017;
  i = i+1;
}
