
#include <iostream>
#include "LocalBranch.hpp"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <queue>

int main()
{
	std::string Folder = "E:\\MipSolver\\collectionunzip";
	std::string filename = Folder + "\\a1c1s1.mps";
	auto scip = ReadFile::readmodel(filename);
	auto first_sol = LocalBranching::getSolution(scip);

	std::queue<std::unordered_map<std::string,double>> sols_q;
	sols_q.push(first_sol);

	int count = 0;
	std::unordered_map<std::string,double> feasible_sol;
	while(!sols_q.empty())
	{

		auto sol = sols_q.front();
		sols_q.pop();
		auto lb_solution = LocalBranching::LocalBranchingFun(scip,sol,20);
		if (lb_solution.first.size() == 1)
		{
			feasible_sol = lb_solution.second;
			break;
		}
		if (lb_solution.first.size() != 0)
		{
			sols_q.push(lb_solution.first);
		}
		if (lb_solution.second.size() != 0)
		{
			sols_q.push(lb_solution.second);
		}
		count++;

		if (count > 1000)
		{
			break;
		}
	}
	std::cout<<"local branching is end "<<std::endl;
	if (feasible_sol.size() > 0)
	{
		std::cout << "find a feasible solution " << std::endl;
	}
	else
	{
		std::cout << "can not  find a feasible solution " << std::endl;
	}

	auto obj = LocalBranching::computeObj(scip, feasible_sol);
	std::cout << "feasible obj is :  " << obj << std::endl;
	return 0;
}
