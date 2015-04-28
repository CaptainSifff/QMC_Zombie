struct Point2
{
double gp;//the Grid-Point. It is the n'th eigenvalue of the matrix
double weight;//this is the first component of the vector
};

static const Point2 grid[] = {//128 points and its weights
{ 2.38748889209e-05  ,  13332.0799669 },
{ 7.16175037393e-05  ,  1481.19399503 },
{ 0.000119338625188  ,  533.122963449 },
{ 0.000167023909988  ,  271.91952327 },
{ 0.000214658993177  ,  164.428172429 },
{ 0.00026222947917  ,  110.015986762 },
{ 0.000309720932555  ,  78.7206753365 },
{ 0.000357118868524  ,  59.0854981328 },
{ 0.000404408742865  ,  45.962749767 },
{ 0.0004515759414  ,  36.7610905374 },
{ 0.000498605768778  ,  30.060670891 },
{ 0.000545483436485  ,  25.0306195045 },
{ 0.000592194049947  ,  21.1584481545 },
{ 0.000638722594565  ,  18.1141258774 },
{ 0.000685053920516  ,  15.6773352003 },
{ 0.000731172726101  ,  13.6964459184 },
{ 0.000777063539429  ,  12.0643322368 },
{ 0.000822710698137  ,  10.7035700191 },
{ 0.000868098326846  ,  9.55707517215 },
{ 0.000913210311952  ,  8.58201187931 },
{ 0.00095803027331  ,  7.74572766511 },
{ 0.00100254153225  ,  7.02297966439 },
{ 0.00104672707527  ,  6.39400374859 },
{ 0.00109056951259  ,  5.84314598644 },
{ 0.00113405103056  ,  5.35787673531 },
{ 0.0011771533367  ,  4.92806976621 },
{ 0.00121985759584  ,  4.54546796464 },
{ 0.00126214435532  ,  4.20328232791 },
{ 0.00130399345689  ,  3.89588748527 },
{ 0.00134538393194  ,  3.61858797296 },
{ 0.001386293876  ,  3.36743695113 },
{ 0.00142670029679  ,  3.13909416858 },
{ 0.00146657892856  ,  2.93071353958 },
{ 0.00150590400242  ,  2.73985319137 },
{ 0.00154464795876  ,  2.56440259799 },
{ 0.00158278108209  ,  2.40252264494 },
{ 0.00162027103037  ,  2.25259530643 },
{ 0.00165708221816  ,  2.11318014582 },
{ 0.00169317499427  ,  1.98297512484 },
{ 0.00172850452679  ,  1.86077930477 },
{ 0.00176301927125  ,  1.74545521046 },
{ 0.00179665886209  ,  1.63589003148 },
{ 0.00182935130274  ,  1.53096197003 },
{ 0.00186100973577  ,  1.42954808504 },
{ 0.00189153101691  ,  1.33072103404 },
{ 0.0019208055358  ,  1.23460796462 },
{ 0.00194876677639  ,  1.14479682754 },
{ 0.00197552265602  ,  1.07124992513 },
{ 0.00200150447436  ,  1.02543857958 },
{ 0.00202735980432  ,  1.0062647343 },
{ 0.00205359832017  ,  1.00107789584 },
{ 0.00208045519613  ,  1.00013513814 },
{ 0.00210801235258  ,  1.00001281315 },
{ 0.00213630795136  ,  1.00000094194 },
{ 0.00216537337487  ,  1.00000005457 },
{ 0.00219524059435  ,  1.00000000252 },
{ 0.00222594326003  ,  1.00000000009 },
{ 0.00225751692329  ,  1.0 },
{ 0.00228999918118  ,  1.0 },
{ 0.00232342982616  ,  1.0 },
{ 0.00235785100877  ,  1.0 },
{ 0.00239330741492  ,  1.0 },
{ 0.00242984645942  ,  1.0 },
{ 0.00246751849755  ,  1.0 },
{ 0.00250637705657  ,  1.0 },
{ 0.00254647908947  ,  1.0 },
{ 0.00258788525353  ,  1.0 },
{ 0.00263066021639  ,  1.0 },
{ 0.00267487299314  ,  1.0 },
{ 0.00272059731781  ,  1.0 },
{ 0.00276791205377  ,  1.0 },
{ 0.00281690164764  ,  1.0 },
{ 0.00286765663229  ,  1.0 },
{ 0.00292027418517  ,  1.0 },
{ 0.00297485874938  ,  1.0 },
{ 0.00303152272556  ,  1.0 },
{ 0.0030903872445  ,  1.0 },
{ 0.00315158303152  ,  1.0 },
{ 0.00321525137559  ,  1.0 },
{ 0.00328154521839  ,  1.0 },
{ 0.00335063038088  ,  1.0 },
{ 0.00342268694821  ,  1.0 },
{ 0.00349791083718  ,  1.0 },
{ 0.0035765155751  ,  1.0 },
{ 0.00365873432395  ,  1.0 },
{ 0.0037448221904  ,  1.0 },
{ 0.00383505886968  ,  1.0 },
{ 0.00392975168128  ,  1.0 },
{ 0.00402923906562  ,  1.0 },
{ 0.00413389462576  ,  1.0 },
{ 0.00424413181578  ,  1.0 },
{ 0.00436040939978  ,  1.0 },
{ 0.00448323783357  ,  1.0 },
{ 0.00461318675629  ,  1.0 },
{ 0.00475089382364  ,  1.0 },
{ 0.00489707517206  ,  1.0 },
{ 0.00505253787593  ,  1.0 },
{ 0.00521819485547  ,  1.0 },
{ 0.00539508281667  ,  1.0 },
{ 0.00558438396814  ,  1.0 },
{ 0.00578745247607  ,  1.0 },
{ 0.00600584690913  ,  1.0 },
{ 0.00624137031733  ,  1.0 },
{ 0.0064961201262  ,  1.0 },
{ 0.00677255076987  ,  1.0 },
{ 0.00707355302631  ,  1.0 },
{ 0.00740255549265  ,  1.0 },
{ 0.00776365576058  ,  1.0 },
{ 0.00816179195343  ,  1.0 },
{ 0.00860296989686  ,  1.0 },
{ 0.00909456817668  ,  1.0 },
{ 0.00964575412678  ,  1.0 },
{ 0.0102680608446  ,  1.0 },
{ 0.0109762029719  ,  1.0 },
{ 0.0117892550438  ,  1.0 },
{ 0.0127323954474  ,  1.0 },
{ 0.0138395602689  ,  1.0 },
{ 0.0151576136278  ,  1.0 },
{ 0.0167531519044  ,  1.0 },
{ 0.018724110952  ,  1.0 },
{ 0.0212206590789  ,  1.0 },
{ 0.0244853758603  ,  1.0 },
{ 0.0289372623803  ,  1.0 },
{ 0.0353677651315  ,  1.0 },
{ 0.0454728408834  ,  1.0 },
{ 0.0636619772368  ,  1.0 },
{ 0.106103295395  ,  1.0 },
{ 0.318309886184  ,  1.0 }
};