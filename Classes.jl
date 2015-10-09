module Classes

type Genotype
  sequence::Vector{Int64}
  num_resistant::Int
  num_drivers::Int
  death_prob_vec::Array{Float64,2}
  growth_prob_vec::Array{Float64,2}
  migration_prob_vec::Array{Float64,2}
  num_cells::Int
  prev_generations
  index
end

# struct Cell {
#   short unsigned int lesion ;
#   short int x,y,z ;
#   unsigned int gen ;
# };

type Location
  loc::Array{Float64,3}
#   x = loc[0] Fix this
#   y = loc[1]
#   z = loc[2]
end


type Cell
  lesion
  location::Location
  genotype::Genotype
end

type Lesion
  int wx
  #Literally rename everything
  r
  rold
  rinit
  rad
  rad0
  n::int  #rename this
  n0::int #rename this
  closest::Vector{Int64}
  nl::int  #rename this
  maxdisp::Float64
end

end
