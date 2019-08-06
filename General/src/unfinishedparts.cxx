
/// TODO

   //   Eigen::HouseHolderQR<MatrixXf> qr(p);
   //Eigen::MatrixXf q = qr.householderQ();
   //std::cout << "Q Matrix :\n" << q << std::endl << std::endl ;
  
   int x = 5;

   //# pragma omp parallel num_threads(thread_count) private(x)
     //{
   // int my_rank = omp_get_thread_num();
   // std::cout << "Thread %d ";
   // x = 2 * my_rank + 2;
   // }
  
 
  /*
  for (int phase = 0; phase < n; phase++) {
    if (phase % 2 == 0) {
# pragma omp parallel for num_threads(thread_count) default(none) shared(a, n) private(i, tmp)
      for (int i = 1; i < n; i += 2) {
	if (a[i-1] > a[i]) {
	  tmp = a[i-1];
	  a[i-1] = a[i];
	  a[i] = tmp;
	}
      }
    } else {
# pragma omp parallel for num_threads(thread_count) default(none) shared(a, n) private(i, tmp)
      for (int i = 1; i < n-1; i += 2) {
	if (a[i] > a[i+1]) {
	  tmp = a[i+1];
	  a[i+1] = a[i];
	  a[i] = tmp;
	}
      }
    } 
  }
  */
  // AA.addup();
  // test_debug();
