
#Method 1: simulation | sample
method_1 <- function(false_positive_rate, population, test_num) {
  results <- rep(0,population-1)
  result_list <- rep(0,test_num)
  for (i in 1:test_num)
  {
    results[sample(population,10)] <- 1
    result_list[i] <- mean(results)
  }
  return(result_list)
}

#Method 2: simulation | rbinomial
method_2 <- function(false_positive_rate, population, test_num) {
  results <- rep(0,population-1)
  result_list <- rep(0,test_num)
  for (i in 1:test_num)
  {
    results <- results | rbinom(population-1,1,prob = false_positive_rate)
    result_list[i] <- mean(results)
  }
  return (result_list)
}

#Method 3: calculation | expected values
method_3 <- function(false_positive_rate, test_num) {
  return(1-(1-false_positive_rate)^c(1:test_num))
}
result_list_1 <- method_1(0.0001, 100000, 10000)
result_list_2 <- method_2(0.0001, 100000, 10000)
result_list_3 <- method_3(0.0001, 10000)

subpoints <- floor(seq(0,100000,length.out = 150))[-1]
plot(subpoints,100*result_list_1[subpoints], main = "Figure. Percentage of Total Population with a False-Positive Test Result", 
     pch = 20, ylim = c(0,70), xlab = "No. of Independent Tests", ylab = "Percentage of Total Population with a False-Positive Test Result")
points(subpoints,100*result_list_2[subpoints], pch = "O", col = "red")
points(subpoints,100*result_list_3[subpoints], col = "green", ylim = c(0,70))
legend("bottomright",c("Method 1","Method 2", "Method 3"), col = c("black","red","green"), pch = c(20,1,0), cex = 0.9)

cor(result_list_1,result_list_3)



      
