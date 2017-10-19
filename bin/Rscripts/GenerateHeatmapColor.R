GenerateHeatmapColor <- function(colFlag, degree){
  my.col = greenred(10)
  if(colFlag == 'bluered'){
    my.col = bluered(10)
  }
  
  col.num = rep(10,9)
  if(degree == 0){
    #nothing
  }else if(degree == 1){
    col.num = c(15,13,12,10, 8,10,12,13,15)
  }else if(degree == 2){
    col.num = c(20,15,12, 8, 4, 8,12,15,20)
  }else if(degree == 3){
    col.num = c(30,15,10, 6, 3, 6,10,15,30)
  }else if(degree == 4){
    col.num = c(40,12, 8, 5, 2, 5, 8,12,40)
  }else if(degree >= 5){
    col.num = c(50,10, 6, 3, 1, 3, 6,10,50)
  }else if(degree == -1){
    col.num = c( 8,10,12,13,15,13,12,10, 8)
  }else if(degree == -2){
    col.num = c( 4, 8,12,15,20,15,12, 8, 4)
  }else if(degree == -3){
    col.num = c( 3, 6,10,15,30,15,10, 6, 3)
  }else if(degree == -4){
    col.num = c( 2, 5, 8,12,40,12, 8, 5, 2)
  }else if(degree == -5){
    col.num = c( 1, 3, 6,10,50,10, 6, 3, 1)
  }else{
    # nothing
  }
  
  color.out = NULL
  for (i in 1:length(my.col)){
    if (i == 1){
      col.tmp = colorpanel(col.num[i],low=my.col[i],high=my.col[i+1]);
      color.out = c(color.out,col.tmp)
    }else if (i == 2){
      col.tmp = colorpanel(col.num[i],low=my.col[i],high=my.col[i+1]);
      color.out = c(color.out,col.tmp)
    }else if (i == 3){
      col.tmp = colorpanel(col.num[i],low=my.col[i],high=my.col[i+1]);
      color.out = c(color.out,col.tmp)
    }else if(i == 4){
      col.tmp = colorpanel(col.num[i],low=my.col[i],high=my.col[i+1]);
      color.out = c(color.out,col.tmp)
    }else if(i == 5){
      col.tmp = colorpanel(col.num[i],low=my.col[i],high=my.col[i+1]);
      color.out = c(color.out,col.tmp)
    }else if(i == 6){
      col.tmp = colorpanel(col.num[i],low=my.col[i],high=my.col[i+1]);
      color.out = c(color.out,col.tmp)
    }else if(i == 7){
      col.tmp = colorpanel(col.num[i],low=my.col[i],high=my.col[i+1]);
      color.out = c(color.out,col.tmp)
    }else if(i == 8){
      col.tmp = colorpanel(col.num[i],low=my.col[i],high=my.col[i+1]);
      color.out = c(color.out,col.tmp)
    }else if(i == 9){
      col.tmp = colorpanel(col.num[i],low=my.col[i],high=my.col[i+1]);
      color.out = c(color.out,col.tmp)
    }else{
      # i = 10
      # doing nothing
    }
  }
  
  return(as.character(color.out))
}

