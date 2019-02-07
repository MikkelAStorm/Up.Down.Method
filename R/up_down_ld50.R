#' @title Analyze Up-Down Lethal Dose 50\% Experiments
#'
#' @description This is an R adaption of the Up-Down method for minimizing experiment subjects in 50\% lethal dose studdies, described by (Brownlee1953). As of this version the code is applied with a commulative likelihood function based on the normal distribution, with exact stepsizes (Brownlee1953) or average stepsizes (flexible or constant delta) as described in (Dixon1980), but not limmeted by the 6-9 expermints as found in Table 7.
#'
#'  Brownlee, K. A., et al. "The Up-and-Down Method with Small Samples." Journal of the American Statistical Association, vol. 48, no. 262, 1953, pp. 262–277., www.jstor.org/stable/2281287.
#'
#'  Dixon, W. J. "Efficient Analysis of Experimental Observations" Annual Review of Pharmacology and Toxicology 1980 20:1, 441-462
#'
#' @param end_point The value of the last experiment in the series.
#' @param steps A vector containing the steps that are used in the experimental series.
#' @param response_pattern A string of X and O, response and no response respectively, representing the chain of results in the up-down experiment chain.
#' @param delta_type Deciding which step type to use. flexible refers to average delta between the used steps, constant is the average delta between all steps in the steps vector, exact is using the exact delta for each step.
#'
#' @return NULL
#'
#' @examples up_down_ld50(end_point = 0.602,
#'     steps = c( 0, 0.301, 0.602, 0.903, 1.204 ),
#'     response_pattern = "oxxoxo",
#'     delta_type = "flexible")
#'
#' @export up_down_ld50 res_to_bol res_bol end_from_response index_from_res flexible_d end_from_response nll.cnorm nll.setD all.character.same
#### working function ####

up_down_ld50 = function(end_point, steps, response_pattern,
                        delta_type = c("exact", "flexible", "constant")){

  # only use first value of the delta type
  delta_type = delta_type[1]
  # Correct response pattern for case
  response_pattern = toupper(response_pattern)

  # check if end_point is one of the steps
  if ( ! end_point %in% steps  ){
    return("End_point not in steps")
  } else {
    # save index of endpoint
    i = match(end_point, steps)
  }

  # check if the response pattern is within the steps
  index = index_from_res( end_index = i, res = response_pattern)

  # check if the response pattern stays within the step ratio
  if (max(index) > length(steps) ){
    return("Response pattern starts above step ratio")
  }
  if (min(index) < 1 ){
    return("Response pattern starts below step ratio")
  }

  # calculate estimated LD50
  # If all characters are the same it is only possible to get a limit
  # and not a precise estimate of the mean.
  if ( all.character.same(response_pattern) ){
    formula_used = F
    char = strsplit(response_pattern, "")[[1]][1]
    # If all X (response) the real mean must be below the end point
    if (char == "X"){
      mu = paste0( "<" , as.character( end_point ) )
    } else{
      # If all O (no response) the mean must be above the end point.
      mu = paste0( ">" , as.character( end_point ) )
    }
    names(mu) = "mean"
    return(mu)
    # use exact delta values
  } else if( delta_type[1] == "exact"){
    # calculate expected mean based on actual steps and negative log likelihood
    # save that the formula is not used
    formula_used = F
    # use average of used delta values as estimate of range where the mean is going to be
    d = flexible_d( end_index = i, res = response_pattern, steps = steps )
    # optimize for mu
    mu = ( optimize(nll.cnorm,
                        end_point + d * c(-2,2),
                        x = steps[index],
                        response =  response_pattern ,
                        sd = d)$minimum )


  } else {
    # calculate based on formula (Dixon1980)
    # Save that the formual is used
    formula_used = T
    # get k from negative log likelihood function, based on the response_pattern:
    k = ( optimize(nll.setD,
                   end_from_response(response_pattern) + c(-10,10),
                   x0 = end_from_response(response_pattern),
                   response = response_pattern )$minimum -
            end_from_response(response_pattern) )

    # get delta value based on delta type
    if (delta_type[1] == "flexible"){
      d = flexible_d(end_index = i, res = response_pattern, steps = steps)

    } else if (delta_type == "constant") {
      d = mean( steps[2 : length(steps)] - steps[1 : (length(steps) - 1)])

    }
    # calcualte mu based on formula using k and delta.
    mu = end_point + d*k
  }

  # return results
  if(formula_used){
    result = as.numeric( c(mu, k, d) )
    names(result) = c("mean", "k", "delta")
  } else{
    result = as.numeric( c(mu) )
    names(result) = c("mean")
  }
  return( result  )
}

#### preperation functions ####

# get boleen response from single X or O
res_to_bol = function( res ){
  # correct to capital letters
  res = toupper(res)
  # if there is response
  if (res == "X"){
    r = T
    # if there is no response
  } else if (res == "O"){
    r = F
  } else {
    return("Unknown response type, can only be X (response) or O (no response)")
    break
  }
  return(r)
}

# get response boolean vector from XO marcow chain
res_bol = function(res){
  res = rapply( as.list( strsplit( res, split = "" )[[1]]) , res_to_bol)
  return(res)
}

# get indexex of steps used based on the response pattern
index_from_res = function( end_index, res){
  # get bool response vector
  res = res_bol(res)
  # initiate index list
  index = rep(0, length(res))
  # set last value
  index[length(index)] = end_index
  # run tough bool response vector to trace back the steps used
  ind = end_index
  for ( i in (length(res)-1): 1){
    if( res[i] ){ ind = ind + 1} else { ind = ind - 1 }
    index[i] = ind
  }
  return(index)
}

# claculate the flexible delta value, based only on the steps used.
flexible_d = function( end_index, res, steps ){
  # get index of steps used
  index = index_from_res( end_index = end_index, res = res)
  # calculate average step size
  d = mean( steps[ (min(index) + 1) : max(index) ] - steps[ min(index) : (max(index) - 1)])
  return(d)
}

# get end position based on response marcow chain.
# It is possible to change the stepsize (d) or start point.
end_from_response = function(response, d = 1, start = 0){
  # start from start point
  end = start
  # get boolean vector from response marcow chain
  res = rapply( as.list( strsplit( response, split = "" )[[1]]) , res_to_bol)
  # use step size and boolean vector to get end position
  for (b in res){
    if( b ){end = end - d} else { end = end + d}
  }
  return(end)
}



# negative log likelihood for the commulative normal distribution.
# Used to optimize mu for likelihood based on response marcow chain
nll.cnorm = function( mu = 0, x, response, sd = 1 ){
  # get boolean vector
  res = rapply( as.list( strsplit( response, split = "" )[[1]]) , res_to_bol)
  # calculate nll for all responses in marcow chain based on given mu.
  nll = 0
  for ( i in 1:length(res)) {
    nll = nll - pnorm( q = x[i], mean = mu, sd = sd, log.p = T, lower.tail = res[i])
  }
  return(nll)
}


# calculate negative log likelihood with a set step size d.
nll.setD = function( mu, d = 1, x0, response ){
  i = 0
  res = rapply( as.list( strsplit( response, split = "" )[[1]]) , res_to_bol)

  # recreate d values for all tests
  D = rep(0, length(res))
  D[length(D) - i] = x0
  for ( b in res[ (length(res) - 1) : 1 ] ){
    i = i+1
    # go down on response, up if no response
    if(b){
      D[length(D) - i] = D[length(D) - i + 1] + d
    } else{
      D[length(D) - i] = D[length(D) - i + 1] - d
    }
  }

  # calculate nll for the d values
  return( nll.cnorm( mu = mu, x = D, response = response)  )
}

# check if all charracters in a string are the same
all.character.same = function(string){
  # split string
  s = strsplit(string, "")[[1]]
  # compare all characters to the one before them
  if(length(s)>1){
    for (i in 1:(length(s)-1 ) ){
      if (s[i] != s[i+1]){
        return (FALSE)
      }
    }
  }
  return (TRUE)
}



