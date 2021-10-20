
for (ell in 7:28){
  cat(all.equal(result_local[[ell]], result_server[[ell]]))
}
