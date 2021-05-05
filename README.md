Les fichiers de ce bloc sont les implémentations de notre mémoire "Simulating Mixtures WithNegative Weights"

Le fichier "Alternate_mixture_simulation.R" permet de simuler un mélange à coefficient alterné de lois normales dans certains cas et est une amélioration par rapport à l'algorithme de Bignami et De Matteis. L'exemple programme est celui de 10000 de variable du mélange 0.5*N(1,1)+1.7*N(0,1)-0.7*N(0,1)-0.5*N(1,1)

Le fichier "Ziggourat-simulation-V7.R" est l'implémentation de notre méthode ziggourats au cas de mélanges de lois normales. Cela permet de simuler 100000 variables suivants la loi tronque 2*N(2,2) -N(2,1) et peut être rapidement modifié pour générer n'importe quel mélange.
