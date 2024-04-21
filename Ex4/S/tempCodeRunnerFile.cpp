if(IsUniform == true){
          for (int i = 0; i < N + 1; ++i)
            r[i] = (i/static_cast<double>(N))*R;
    }else{
          for(int i = 0; i < N+1 ; ++i)
            r[i] = sqrt(i/static_cast<double>(N))*R ; 
    };