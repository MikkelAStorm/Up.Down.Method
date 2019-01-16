#' @title Analyze Up-Down Lethal Dose 50\% Experiments
#'
#' @description This is an R adaption of the Up-Down method for minimizing experiment subjects in 50\% lethal dose studdies, described by (Brownlee1953). It is developed for Von Frey filament experiments, thus standard steps for these are hardcoded, but with the standard costum sitting it is applicable for all Up-Down-ld50 expreiments. As of version 0.1.0 the code is applied with a commulative likelihood function based on the normal distribution, with exact stepsizes (Brownlee1953) or average stepsizes as described in (Dixon1980), but not limmeted by the 6-9 expermints as found in Table 7. Brownlee, K. A., et al. "The Up-and-Down Method with Small Samples." Journal of the American Statistical Association, vol. 48, no. 262, 1953, pp. 262–277., www.jstor.org/stable/2281287. Dixon, W. J. "Efficient Analysis of Experimental Observations" Annual Review of Pharmacology and Toxicology 1980 20:1, 441-462
#'
#' @param end_point The weight of the last experiment in the series.
#' @param steps A vector containing the steps that are used in the experimental series.
#' @param response_pattern A string of X and O, response and no response respectivly, repressenting the chain of results in the up-down experiment chain.
#' @param ID ID of the experiment chian - used to differentiate the output for multiple tests
#' @param log_steps This flag is to differentiate between if the steps are original values, or log 10 scale)
#' @param delta_type Diciding which step type to use. flexible referes to average delta between the used steps, constant is the average delta between all steps in the steps vector, exact is using the exact delta for each step.
#' @param verbose
#'
#' @return NULL
#'
#' @examples up_down_ld50(end_point = 0.07,
#'     steps = c( 0.008, 0.02, 0.04, 0.07, 0.16, 0.4, 0.6, 1, 2),
#'     log_steps = F,
#'     response_pattern = "xoxoxx", ID = 1,
#'     delta_type = "flexible")
#'
#' @export up_down_ld50 res_to_bol res_bol end_from_response index_from_res flexible_d end_from_response nll.cnorm nll.setD all.character.same
#### working function ####

up_down_ld50 = function(end_point, steps, log_steps = FALSE,
                        response_pattern, ID, verbose = FALSE,
                        delta_type = c("flexible", "constant", "exact")){

  # only use first value in each step
  delta_type = delta_type[1]
  response_pattern = toupper(response_pattern)
  steps_text = paste(steps, collapse=', ')

  # check if end_point is one of the steps
  if ( ! end_point %in% steps  ){
    return("end_point not in steps")
  } else {
    i = match(end_point, steps)
  }

  # change steps to to log scale
  if ( ! log_steps){
    steps = log10(steps)
  }

  # calculate estimated LD50
  if ( all.character.same(response_pattern) ){
    char = strsplit(response_pattern, "")[[1]][1]
    if (char == "X"){
      mu = paste0( "<" , as.character( end_point ) )
    } else{
      mu = paste0( ">" , as.character( end_point ) )
    }
  } else if( delta_type[1] == "exact"){
    if ( ! log_steps ){
      end_point = log10(end_point)
    }
    # calculate expected mean based on actual steps and negative log likelihood
    index = index_from_res( end_index = i, res = response_pattern)
    d = flexible_d( end_index = i, res = response_pattern, steps = steps )
    # optimize for mu
    mu_opt = ( optimize(nll.cnorm,
                        end_point + d * c(-2,2),
                        x = steps[index],
                        response =  response_pattern ,
                        sd = d)$minimum )

    if (log_steps){
      mu = as.character(10^(mu_opt)/10000)
    } else{
      mu = as.character(10^(mu_opt))
    }


  } else {

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

    if (log_steps){
      mu = as.character( 10^( end_point + d*k )/10000 )
    } else {
      mu = as.character( 10^( log10( end_point ) + d*k ) )
    }

  }

  if (verbose){
    if(log_steps){
      method_expl = "using the handle number"
    } else {
      method_expl = "force"
    }

    text = sprintf("Sensitivity to evoked mechanical stimulation was assessed using a series of calibrated von Frey monoﬁlaments (%s). The Dixon up and down method was used to calculate the 50 percent withdrawal threshold [1]. The calculation was made with up_down_ld50 {up.down.ld50} [2] %s and a %s delta-value [3].
    1. Dixon, W. J. 'Efficient Analysis of Experimental Observations' Annual Review of Pharmacology and Toxicology 1980 20:1, 441-462
    2. Reference to program
    3. Reference to our manuscript", steps_text, as.character(method_expl), as.character(delta_type ))

    return (text[1])
  } else {
    return (mu)
  }
}



#### preperation functions ####
# get response boolean vector from XO marcow chain
res_to_bol = function( res ){
  # correct for normal or capital letters
  res = toupper(res)
  # if there is response
  if (res == "X"){
    r = T
    # if there is no response
  } else if (res == "O"){
    r = F
  } else {
    print("Unknown response type, can only be X (response) or O (no response)")
    break
  }
  return(r)
}

res_bol = function(res){
  res = rapply( as.list( strsplit( res, split = "" )[[1]]) , res_to_bol)
  return(res)
}

index_from_res = function( end_index, res){
  res = res_bol(res)
  index = rep(0, length(res))
  index[length(index)] = end_index
  ind = end_index
  for ( i in (length(res)-1): 1){
    if( res[i] ){ ind = ind + 1} else { ind = ind - 1 }
    index[i] = ind
  }
  if( any( index < 1 ) ){
    print("index below 1")
    break
  }
  return(index)
}

flexible_d = function( end_index, res, steps ){
  index = index_from_res( end_index = end_index, res = res)
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

all.character.same = function(string){
  s = strsplit(string, "")
  s = s[[1]]
  for (i in 1:(length(s)-1 ) ){
    if (s[i] != s[i+1]){
      return (FALSE)
    }
  }
  return (TRUE)
}



