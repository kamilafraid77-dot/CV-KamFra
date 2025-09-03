############################################################
# TASK 1: Pathfinding in a Maze
############################################################
# Algorithm description:
# This program loads a maze stored in an RDS file, visualizes it, 
# and checks if there exists a valid path from a starting point 
# (top-left corner) to a target region (defined as the "logo" area).
#
# It uses a breadth-first search (BFS) approach:
# 1. Reads and displays the maze as a raster image.
# 2. Defines a function `pathQ` that:
#    - Validates coordinates inside the maze.
#    - Iteratively explores neighboring cells using a queue.
#    - Keeps track of visited cells to avoid cycles.
#    - Returns TRUE if the target region is reached, otherwise FALSE.
# 3. Tests the algorithm by checking if there is a path from (1,1) 
#    to the defined end region and prints the result.
############################################################

file_path <- file.choose()
d0 <- readRDS(file = file_path)
plot(as.raster(d0))

pathQ <- function(maze, startPoint, endRegion) {
  isValid <- function(x, y) {
    x >= 1 && x <= ncol(maze) && y >= 1 && y <= nrow(maze)
  }
  
  start_x <- startPoint$x
  start_y <- startPoint$y
  end_x <- endRegion$x
  end_y <- endRegion$y
  
  queue <- list(list(x = start_x, y = start_y))
  visited <- matrix(FALSE, nrow = nrow(maze), ncol = ncol(maze))
  
  while (length(queue) > 0) {
    current <- queue[[1]]
    queue <- queue[-1]
    
    if (current$x %in% end_x && current$y %in% end_y) {
      return(TRUE)  
    }
    
    neighbors <- list(
      list(x = current$x - 1, y = current$y),
      list(x = current$x + 1, y = current$y),
      list(x = current$x, y = current$y - 1),
      list(x = current$x, y = current$y + 1)
    )
    
    for (neighbor in neighbors) {
      if (isValid(neighbor$x, neighbor$y) && 
          maze[neighbor$y, neighbor$x] && 
          !visited[neighbor$y, neighbor$x]) {
        queue <- c(queue, list(neighbor))
        visited[neighbor$y, neighbor$x] <- TRUE
      }
    }
  }
  return(FALSE)
}

### End region
logoPosition <- list(x = 387:413, y = 322:348)

### Starting point
startPoint <- list(x = 1, y = 1)

### Run test
result <- pathQ(d0, startPoint = startPoint, endRegion = logoPosition)
print(result)


############################################################
# TASK 2: Pure Strategy Nash Equilibria
############################################################
# Algorithm description:
# This program defines a function `getAllPureStrategyNE` that computes 
# all pure strategy Nash equilibria in a normal-form game.
#
# Steps:
# 1. Takes as input a game represented as a list of payoff matrices 
#    (one per player).
# 2. Iterates over all possible action profiles (combinations of strategies).
# 3. For each profile, checks whether any player has an incentive to deviate:
#    - Compares the player's payoff in the current profile with payoffs 
#      from all possible unilateral deviations.
#    - If a player can improve by deviating, the profile is not an equilibrium.
# 4. If no player benefits from deviating, the profile is stored as a pure 
#    strategy Nash equilibrium.
# 5. Returns a list of all such equilibria.
#
# Example: a simple 2-player coordination game where the function 
# identifies pure strategy Nash equilibria.
############################################################

getAllPureStrategyNE <- function(gra) {
  num_graczy <- length(gra)
  num_akcji <- sapply(gra, function(gracz) dim(gracz)[1])
  
  pure_strategy_ne <- list()
  
  for (i in 1:num_akcji[[1]]) {
    for (j in 1:num_akcji[[2]]) {
      biezacy_profil <- c(i, j)
      rownowaga_nasha <- TRUE 
      
      for (gracz in 1:num_graczy) {
        akcje_innych_graczy <- biezacy_profil[-gracz]
        rozszerz_akcje_innych_graczy <- expand.grid(
          replicate(num_graczy - 1, 1:num_akcji[[gracz]], simplify = FALSE)
        )
        
        for (action in 1:num_akcji[[gracz]]) {
          rozszerz_akcje_innych_graczy[[gracz]] <- action
          inny_profil <- as.numeric(unlist(rozszerz_akcje_innych_graczy))
          
          if (gra[[gracz]][biezacy_profil[gracz]] < gra[[gracz]][inny_profil[gracz]]) {
            rownowaga_nasha <- FALSE
            break
          }
        }
        if (!rownowaga_nasha) break
      }
      
      if (rownowaga_nasha) {
        pure_strategy_ne[[paste(biezacy_profil, collapse = "")]] <- biezacy_profil
      }
    }
  }
  return(pure_strategy_ne)
}

### Example (coordination game for two players)
game <- list(
  "player1" = array(c(1, 0, 0, 1), dim = c(2, 2)),
  "player2" = array(c(1, 0, 0, 1), dim = c(2, 2))
)
getAllPureStrategyNE(game)


############################################################
# TASK 3: Monte Carlo Integration
############################################################
# Algorithm description:
# This program defines a function `integrate3d` that estimates a 2D integral 
# of a given function f(x, y) over a rectangular region using the Monte Carlo method.
#
# Steps:
# 1. Takes as input:
#    - f: the function to integrate,
#    - over: a list with ranges for x and y,
#    - n: the number of random samples.
# 2. Repeatedly generates random points (x, y) uniformly distributed in the region.
# 3. For each point, checks whether a randomly drawn value is below f(x, y). 
#    - This simulates sampling under the surface of the function.
# 4. Counts the fraction of points that fall below the function surface.
# 5. Multiplies this fraction by the area of the region to approximate the integral.
# 6. Returns the estimated integral value.
#
# Examples:
# - Low number of samples (n = 100) → rough approximation.
# - High number of samples (n = 100,000) → more accurate result.
############################################################

integrate3d <- function(f, over, n) {
  x_przedz <- over$x 
  y_przedz <- over$y
  
  ponizej_wyk <- function(x, y) { 
    f_wart <- f(x, y)
    return(runif(1) < f_wart)
  }
  
  punkty_ponizej <- 0
  
  for (i in 1:n) {
    x <- runif(1, x_przedz[1], x_przedz[2])
    y <- runif(1, y_przedz[1], y_przedz[2])
    
    if (ponizej_wyk(x, y)) {
      punkty_ponizej <- punkty_ponizej + 1
    }
  }
  
  total_wart <- (x_przedz[2] - x_przedz[1]) * (y_przedz[2] - y_przedz[1]) 
  integral_wart <- punkty_ponizej / n * total_wart
  return(integral_wart)
}

### Example 1 (low n)
result_low_n <- integrate3d(
  f = function(x, y) {cos(x) * y},
  over = list(x = c(0, pi / 2), y = c(0, 1)),
  n = 10^2
)

### Example 2 (high n)
result_high_n <- integrate3d(
  f = function(x, y) {cos(x) * y},
  over = list(x = c(0, pi / 2), y = c(0, 1)),
  n = 10^5
)

print(result_low_n)
print(result_high_n)
