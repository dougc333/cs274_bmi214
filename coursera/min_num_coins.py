def RecursiveChange(money,Coins):
  if money==0:
    return 0
  MinNumCoins = 10000
  for i in range(0,len(Coins)-1):
    if money >= Coins[i]:
      numcoins=RecursiveChange(money-Coins[i],Coins)
      if numcoins+1 < MinNumCoins:
        MinNumCoins = numcoins+1
  return MinNumCoins
print(RecursiveChange(76,[5,4,1]))