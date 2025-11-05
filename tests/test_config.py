import yaml

config = yaml.safe_load(open('aws_test/configs/phase2_miller_urey_extended.yaml'))
print('Config loaded successfully!')
print(f"Scenario: {config['scenario_name']}")
print(f"Steps: {config['simulation']['max_steps']}")
print(f"Particles: {config['simulation']['n_particles']}")

