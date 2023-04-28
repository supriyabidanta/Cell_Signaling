cd Signaling/Data

curl -o local.h5ad "https://corpora-data-prod.s3.amazonaws.com/0b4a15a7-4e9e-4555-9733-2423e5c66469/local.h5ad?AWSAccessKeyId=ASIATLYQ5N5X4MOSGHJP&Signature=Ck5WNFMDV9nyUcCnsg%2FnbhbYyLI%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEEAaCXVzLXdlc3QtMiJGMEQCIBS7MeTO5WsFuChrzDT93qI0npDo6mtpHw3YlAlqlcjLAiAM%2BrpwLB7R37bbXLiYCGnu858mKKLme0osUGN2Vl3MpCrrAwhJEAEaDDIzMTQyNjg0NjU3NSIMCykux7qQeB7GWO%2B1KsgDytuIlB9p3JdN%2FAJp1Vx9lIqveFdFreKobBqzJM%2FZ%2BiYQiwDeRarl3UctEYdzTrJUHEWF3gkSpPsbjRaaGk8Pw562me4k%2B3u7g2IbKxi3IgUCP8UbZEzLTuOxWSwcWfIqhT2HPnWOsmmWnz4%2BoFaK4h0A5XyQ%2FtElkRxf3Id5oZvmO0w%2FSj4ArtmK%2F4LhFW%2FQ4X1DKDjm1IhuwBfnFZH0dlTvt7KFaUA3crXWDTIGltSBxCy7d3y4bYI25a6unxCeXD16%2BuELGed3DYAwpfCTgwfz6YykMpj4%2FuHqrDG1Yus2%2Fuo8bFh%2BYlR6RDkqAymgoEKtNTc4G%2FM3nNGM7BZn8%2Bkkr0rYR67WjaQehR37GInST3EX210BbkQIOKYu6qlyCt75Luq%2FB2BFihtFtqQ4NyRDgPKwJI5hs%2BGAqj%2FtXwvaRdeSlHS9PnYKzrEmeavFOy8b%2BRZ91Bhnun2VYA1KSZA3ttXphHuOwJrode9%2FctMMCKBjgw8s7lwMfMglOXhDB8k0JiV3Nb3aB7mqltNflmlnrj97fbvkRCndEh6qWxpJfFAHls%2FiZ2jfkP5uTfnRWbl27bymx%2BHVzSpnBB%2BeSeKyT7ZM1vqrMO6%2BqqIGOqYB%2BSTB1lWYWnuzbH6T0asNsx2TgKGs2FO4JSmnzVWuGnr5055IoMltrLN63siiXDaw4xcl2OtTNPECveumVO6PV%2BVPkSkoQVl3WJki5VOoDsRltZVUSI0p5oNrVGbN1oNiZJ4VLjm9HqD2tnh0NHuW2iQRrMBHIShrvmaahaBTjzVrq1MTtrgQ6io57xWfNfHEfv3k1X7MGs69o6fN%2BeyjUy4BIwX0kA%3D%3D&Expires=1683221818"

curl -o local.rds "https://corpora-data-prod.s3.amazonaws.com/0b4a15a7-4e9e-4555-9733-2423e5c66469/local.rds?AWSAccessKeyId=ASIATLYQ5N5XROTTH3U4&Signature=5%2BaDfoNv0b0Nn7K5sK1%2B4FpmNUs%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEEAaCXVzLXdlc3QtMiJHMEUCIQDyL6vkC7hEq8UTLbmNNly9ZEyD4tEhOd72jNi18NVPaQIgAss3zZqZ89mJ7tO7faQjc6yXYHWdJXcTftqrSZUnqjMq6wMISRABGgwyMzE0MjY4NDY1NzUiDMaimD6qU7R8Ts0lfCrIAzlJgEERFFmMDcSOKQUVZG40DMb6rAaGulDK6Zsh4d9ai0PE%2B8aS7HZqRsIkBV9sIqxG3IyvINiGRTdgIyoWlVBdJtZf2ZPxd4754H8X5Fyl9zXJsa3P1I9pP1HhxiAlcQ9DdMhpwzIGDXYQkrxgLtbtNBlaiqUEqMHgrpPAi0ck4649TRWhh62ToEaitHrqPCLRL%2BjObEzIG6RVxAS1ZtrlU6hrlD%2B8wcf6t7%2FKPV1MeXicYlavuI%2FGLzDuGXO9e%2F5qnf6t0ITI9CpQHUb3gmIJj6ny3AkNgH7Js9fdXnjYTHXzXSQOE1EeKM%2FnZEPdt7qQunw5lHO83o0eTJDMw7vV%2BpnJzA9EKb9Prbykj8bE9M5%2FIFL6iWdx6BWsa5cmtqJ2n%2BcWXyeP75MTUyCd%2BwCVnLZrOBP9JCyG7bp1HerfnK42%2FoQbLyCqJGpIwfQIrmzz6htKeaM7VDjs4wpPVdR%2BnbFFqBi5m50pLbYi2jBcmIUzEj2YxHd8Zx0a6uV%2Fb2n6zfPVoY54EnukDjdsL6ulWLC0CVMYe1zDBdRHLj1EZYgio7uEscBKgWh1txwmrEmnwGD373OdqcU%2BR3dKnOFd%2Bvy%2B1p62HDCMvqqiBjqlAfVQShsXg4K8bJ9x60ksOWHuf1QasDlJX7WIjvuysLUpUJDaZzZSpy%2FEwipRKQr63IeHmrMcELm8AEpBJ4n1tAuRgQ5zOxsfleyLfKfq7tLXXlQ3ESGPprUAPBiQmWKjv1hOQjsXZ%2ByIoQm1JzTVUx0h4Q8soxRSQlt1lhOink0wPbbu1kbFl3OfU1NGIMfjZGfQTjXtiJyXByg1mkpvogao2Zvp8w%3D%3D&Expires=1683221839"


