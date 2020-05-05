*---------------------------------------------------------------------------------------------------------------------------------------#
*chebyshev recursion or chebyshve polynomial generation
$macro chebyrecgentst(i_basis,i_state,grid_state)   tst(i_basis,i_state); \
                                                    loop(i_basis, \
                                                    tst(i_basis,i_state)$(ord(i_basis) eq 1) = 1; \
                                                    tst(i_basis,i_state)$(ord(i_basis) eq 2) = grid_state; \
                                                    tst(i_basis,i_state)$(ord(i_basis) ge 3) = 2*grid_state*tst(i_basis-1,i_state) - tst(i_basis-2,i_state); \
                                                    ); 

$macro fdchebyrecgentst(i_basis,i_state,grid_state)   tst2(i_basis,i_state); \
                                                    loop(i_basis, \
                                                        tst(i_basis,i_state)$(ord(i_basis) eq 1) = 1; \
                                                        tst(i_basis,i_state)$(ord(i_basis) eq 2) = grid_state; \
                                                        tst(i_basis,i_state)$(ord(i_basis) ge 3) = 2*grid_state*tst(i_basis-1,i_state) - tst(i_basis-2,i_state); \
                                                        ); \
                                                    loop(i_basis, \
                                                        tst2(i_basis,i_state)$(ord(i_basis) eq 1) = 0; \
                                                        tst2(i_basis,i_state)$(ord(i_basis) eq 2) = 1; \
                                                        tst2(i_basis,i_state)$(ord(i_basis) ge 3) = ((ord(i_basis)-1)/2)*(1/(1-sqr(grid_state)))*(tst(i_basis-1,i_state) - tst(i_basis+1,i_state) ); \
                                                        tst2(i_basis,i_state)$(ord(i_basis) eq card(i_basis)) = (ord(i_basis)-1)*sin((ord(i_basis)-1)*arccos(grid_state))/(sqrt(1 - (sqr(grid_state)))); \
                                                         ); \

$macro chebyrecgen(i_basis,grid_state)  cos((ord(i_basis)-1)*arccos(grid_state))

*chebyshev recursion or chebyshve polynomial generation for first derivative
$macro fdchebyrecgen(i_basis,grid_state)  (ord(i_basis)-1)*sin((ord(i_basis)-1)*arccos(grid_state))/(sqrt(1 - (sqr(grid_state))))

*chebyshev recursion or chebyshve polynomial generation for second derivative
*$macro sdchebyrecgen(i_basis,grid_state)    (ord(i_basis)-1)*(((ord(i_basis)-1)*cos((ord(i_basis)-1)*arccos(grid_state)))/(-1 + sqr(grid_state)) \
*                                            + (grid_state*sin((ord(i_basis)-1)*arccos(grid_state)))/rpower((1 - sqr(grid_state)),(3/2)))

*$macro sdchebyrecgen(i_basis,grid_state)    sqr(ord(i_basis)-1)*chebyrecgen(i_basis,grid_state)/(-1+sqr(grid_state))  \
*                                            +  grid_state*fdchebyrecgen(i_basis,grid_state)/(1-sqr(grid_state))

$macro sdchebyrecgen(i_basis,grid_state)   (ord(i_basis)-1)*( ((ord(i_basis)-1)*cos((ord(i_basis)-1)*arccos(grid_state)))/(-1+grid_state*grid_state) \
                                            + (grid_state*sin((ord(i_basis)-1)*arccos(grid_state)))/(vcpower((1-grid_state*grid_state),(3/2))) ) 

*chebyshev collocation nodes
$macro chebycolnod(i_state)     cos( pi*((card(i_state) - ord(i_state)) + 0.5)/card(i_state))
*normalized state grid i.e. state grid on [-1,1]
$macro normstategri(grid_state,min_state,max_state)     2*(grid_state - min_state)/(max_state - min_state) - 1
*using cheby collocation nodes to generate (shifted) state grid i.e. grid on [a,b]
$macro shifstategri(colnodes,min_state,max_state) (min_state + max_state)/2 + ((max_state - min_state)/2)*colnodes
*---------------------------------------------------------------------------------------------------------------------------------------#