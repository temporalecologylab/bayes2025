//
// This (mathematically incoherent) Stan program defines a simple ice cream 
// recipe, for geeks and people with too much time on their hands
// (see normal recipe below)
//
//

// The input data is a vector 'happiness' (on a scale from 1 to 10?) 
// of length 'N' (=number of people in bayes meeting)
data {
  int<lower=0> N;
  vector[N] happiness;
}

transformed data {
  // ingredients
  int blackberries = 4; //cups
  real lemonjuice = 0.5; // teaspoons
  real sugar = 0.5; //cups
  int whippingcream = 2; //cups (33% milk fat...or something)
  int sweetenedcondensedmilk = 14; // ounces (I guessed, so...draw from N(14,1)?)
  int vanillaextract = 2; // teaspoons
  
  //supplies
  int blender = 1;
  int sieve = 1;
  int stove = 1;
  int saucepan = 1; // medium sized
  int freezer = 1;
  
  // first steps
  real puree = blackberries*blender; // blend until you get a "smooth, juicy consistency"
  real blackberryjuice = puree + sieve; // press puree through a fine-mesh sieve and discard pulp--should leave ~ 2 cups liquid 
  real syrup = ((blackberryjuice + lemonjuice + sugar + saucepan)/stove) *10; // simmer + stir over medium-high heat until reduce by half and thickened; then transfer to a bowl and refrigerate until cold
  
  real whippedcream = whippingcream*blender; // briefly until foamy?
  real mixedcream = (whippedcream + sweetenedcondensedmilk + vanillaextract)*blender; // until thickened
}

parameters {
  real blackberry_tastiness; // effect of blackberry quality on happiness?
  real cream_tastiness; // effect of whipped cream on happiness?
  real<lower=0> sigma; // variation in how happy people are to be in bayes meeting, generally speaking?
}

model {
  for (n in 1:N){
    target += log_mix(0.5, normal_lpdf(happiness[n] | blackberry_tastiness*syrup, sigma),
                            normal_lpdf(happiness[n] | cream_tastiness*mixedcream, sigma));
  }
  // no clue what this means mathematically but...add two cups of mixed cream to the puree and "mix" until well combined!
  // then gently fold in the remaining cream to create swirls
  
  target += 6*freezer; // hours in the freezer
}

// syntactically correct, but will it compile? who knows!
// (this is an exercise left for the reader)


// Okay now the version that people actually care about:

// Ingredients //

// 24 ounces (680 g) fresh blackberries, rinsed and patted dry, about 4 cups
// 0.5 teaspoon lemon juice
// 0.5 cup (100 g) sugar
// 2 cups (480 ml) heavy cream
// 14 ounces (397 g) sweetened condensed milk
// 2 teaspoons vanilla extract


// Instructions //

// 1
// Add the blackberries to a blender and pulse until fully pureed. You want a 
// smooth, juicy consistency. Press the puree through a fine-mesh sieve and 
// discard the pulp. You should have about 2 cups of liquid.

// 2
// Transfer the liquid to a medium saucepan, add in the lemon juice and sugar. 
// Heat over medium-high heat until reduced by half and thickened, about 5 minutes; 
// stir occasionally. Transfer to a bowl and refrigerate until cold, about 2 hours.

// 3
// Add the heavy cream to a large bowl. Using and electric mixer beat until stiff 
// peaks form. Turn the speed to low and slowly pour in the condensed milk. Add 
// the vanilla extract. Turn the speed back up to high and beat until thickened.

// 4
// Transfer 2 cups of the whipped cream to a separate bowl. Add in the cooled 
// puree and mix until well combined. Gently fold into the remaining whipped 
// cream to create swirls. Don't over mix or you will lose the swirls.


// 5
// Transfer to a freezer-safe container and freeze for at least 6 hours.
