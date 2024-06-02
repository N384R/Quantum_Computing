from spinor import Spinor

class SlaterDeterminant:
    def __init__(self, K):
        self.K = K
        self.orbital = self.get_orbital(K)
    
    def get_orbital(self, K):
        orbital = {}
        for i in range(K):
            orbital[f'X{i}'] = Spinor('1')
        return orbital
    
    def __repr__(self):
        return f'{self.orbital}'
    

if __name__ == "__main__":
    slater = SlaterDeterminant(6)
    print(slater)
    print(type(slater.orbital['X1']))
