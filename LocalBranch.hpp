#pragma once
#include <iostream>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <map>
#include <set>
#include <iostream>
#include <set>
#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <random>
#include <unordered_set>

namespace ReadFile
{
	SCIP* readmodel(std::string filename)
	{
        SCIP* scip = NULL;
        SCIPcreate(&scip);
        SCIPincludeDefaultPlugins(scip);
        // 设置一些求解参数
        SCIPsetIntParam(scip, "presolving/maxrounds", 0);
        SCIPsetIntParam(scip, "presolving/qpkktref/maxrounds", 0);
        SCIPsetBoolParam(scip, "constraints/linear/upgrade/indicator", false);
        //将string类型的filename转换成const char*类型
        auto filename1 = filename.c_str();
        // 读入模型的mps文件
        SCIPreadProb(scip, filename1, NULL);
        int nrows = SCIPgetNConss(scip);
        int ncols = SCIPgetNVars(scip);
        std::cout<<"load model success !"<<"the num of cons is "<<nrows<<"  the num of vars is "<<ncols<<std::endl;
        //SCIPsolve(scip);
        return scip;
	}

    double random()
    {
        // 生成随机数
        double r = (double)rand() / RAND_MAX;
        return r;
    }

    int randomInt(int a,int b)
	{
        std::mt19937 gen(1024); // 使用梅森旋转算法初始化一个随机数生成器
        std::uniform_int_distribution<> dis(a, b); //[a,b]

        // 生成并输出随机数
        return dis(gen);
	}

    int randomDouble(int a, int b)
    {
        std::mt19937 gen(1024); // 使用梅森旋转算法初始化一个随机数生成器
        std::uniform_real_distribution<> dis(a, b);
        // 生成并输出随机数
        return dis(gen);
    }

}

namespace LocalBranching
{
	std::unordered_map<std::string, double> getSolution(SCIP* scip)
	{
        std::unordered_map<std::string, double> tmpsolution;
        SCIP_VAR** vars_ori = SCIPgetVars(scip);
        int ncols = SCIPgetNVars(scip);
        for (int i = 0; i < ncols; i++)
		{
			SCIP_VAR* var = vars_ori[i];
			std::string varname = SCIPvarGetName(var);
            auto var_type = SCIPvarGetType(var);
			if (var_type == SCIP_VARTYPE_BINARY)
			{
                auto r = ReadFile::random();
                if (r > 0.5)
                {
	                tmpsolution[varname] = 1;
				}
				else
				{
					tmpsolution[varname] = 0;
                }
			}
			else if(var_type == SCIP_VARTYPE_INTEGER)
			{
				auto lb = SCIPvarGetLbLocal(var);
                auto ub = SCIPvarGetUbLocal(var);
				if (ub - lb > 1000)
				{
                    tmpsolution[varname] = lb;
				}
				else
				{
                    auto value = ReadFile::randomInt(lb, ub);
                    tmpsolution[varname] = value;
				}
			}
			else
			{
				auto lb = SCIPvarGetLbLocal(var);
				auto ub = SCIPvarGetUbLocal(var);
				if (ub - lb > 1000)
				{
                    tmpsolution[varname] = lb;
				}
				else
				{
                    auto value = ReadFile::randomDouble(lb, ub);
                    tmpsolution[varname] = value;
				}
			}
		}
        
        return tmpsolution;
	}

    double computeObj(SCIP* scip,const std::unordered_map<std::string,double>& sol)
	{
        double obj = 0.0;
        int nvars = SCIPgetNVars(scip);
        SCIP_VAR** vars = SCIPgetVars(scip);

        for (int i = 0; i < nvars; i++)
        {
            std::string var_name = SCIPvarGetName(vars[i]);
            auto it = sol.find(var_name);

            if (it != sol.end())
            {
                double var_value = it->second;
                double obj_coef = SCIPvarGetObj(vars[i]);
                obj += obj_coef * var_value;
            }
        }

        return obj;
	}

    
    std::pair<std::unordered_map<std::string, double>, std::unordered_map<std::string, double>> LocalBranchingFun(SCIP* _scip,const std::unordered_map<std::string, double>& tmpsolution,double neighbor)
	{
        int _ncols = SCIPgetNVars(_scip);
        int _nrows = SCIPgetNConss(_scip);
        std::cout << "the num of vars of model is " << _ncols << "  the num of cons of model is " << _nrows << std::endl;

        std::unordered_set<std::string> binary_fixed;
        std::unordered_set<std::string> binary_free;
        std::unordered_set<std::string> redundant_var;
        std::unordered_map<std::string, SCIP_VAR*> varmap;
        SCIP* scip = NULL;
        SCIPcreate(&scip);
        // 使用 SCIPcopy() 来复制原始问题
        SCIP_Bool valid;
        SCIPcopyOrig(_scip, scip, NULL, NULL, "_copy", TRUE, FALSE, TRUE, &valid);
        SCIPsetMessagehdlrQuiet(scip, true);
        int ncols = SCIPgetNVars(scip);
        int nrows = SCIPgetNConss(scip);
        std::cout<<"the num of vars of copymodel is "<<ncols<<"  the num of cons of copymodel is "<<nrows<<std::endl;

        // 创建一个新的解
        SCIP_SOL* sol;
        SCIPcreateSol(scip, &sol, NULL);

        SCIP_VAR** vars_ori = SCIPgetVars(scip);
        for (int i = 0; i < ncols; ++i)
        {
			SCIP_VAR* var = vars_ori[i];
			std::string varname = SCIPvarGetName(var);
            //std::cout << varname << std::endl;
			varmap[varname] = var;

            //change obj
            SCIP_Real newcoef = 0;
            SCIPchgVarObj(scip,var,newcoef);
        }
        // 为每个变量设置值
        for (const auto& item : tmpsolution)
        {
            SCIPsetSolVal(scip, sol, varmap.at(item.first), item.second);
            auto var_type = SCIPvarGetType(varmap.at(item.first));
            if (var_type == SCIP_VARTYPE_BINARY)
            {
	            if (item.second == 1)
	            {
		            binary_fixed.insert(item.first);
	            }
	            else
	            {
		            binary_free.insert(item.first);
	            }
            }
        }
        
        // 检查每个约束

        
        SCIP_CONS** conss = SCIPgetConss(scip);
        for (int i = 0; i < nrows; ++i)
        {
            SCIP_CONS* consdata = conss[i];
            SCIP_RESULT result;
        	SCIPcheckCons(scip,consdata,sol,TRUE,TRUE,TRUE,&result);
			if (result == SCIP_INFEASIBLE)
			{
                double tmpvalue = 0;
                double delta = 0;
                // 获得约束的系数和相应的变量
                int nnonz = SCIPgetNVarsLinear(scip, consdata);
                SCIP_VAR** vars = SCIPgetVarsLinear(scip, consdata);
                SCIP_Real* coefs = SCIPgetValsLinear(scip, consdata);
                // 遍历一条约束中的非零元素
                for (int j = 0; j < nnonz; ++j)
                {
                    auto temp_name = SCIPvarGetName(vars[j]);
                    //std::cout<<"coef is "<<coefs[j]<<"  x is  "<< tmpsolution.at(temp_name) <<std::endl;
                    tmpvalue += coefs[j] * tmpsolution.at(temp_name);
                }
                unsigned success;
                auto Lb = SCIPconsGetLhs(scip, consdata, &success);
                auto Ub = SCIPconsGetRhs(scip, consdata, &success);
                std::string consname = SCIPconsGetName(consdata);
                if (Lb==Ub)
                {
                    delta = Ub - tmpvalue;
                    //std::cout << Lb << "  " << tmpvalue << "  " << delta << std::endl;
                }
                else if (Lb > -SCIPinfinity(scip))
                {
	                delta = Lb - tmpvalue;
                    //std::cout << Lb << "  " << tmpvalue << "  " <<delta<< std::endl;
                }
				else
				{
	                delta = Ub - tmpvalue;
                    //std::cout << Ub << "  " << tmpvalue << "  " << delta << std::endl;
				}

                //add variable
                SCIP_VAR* y;
                std::string yname = "y_" + consname;
                SCIPcreateVarBasic(scip, &y,yname.c_str() , 0, 1, 1.0, SCIP_VARTYPE_BINARY);
                SCIPaddVar(scip,y);
                SCIPaddCoefLinear(scip, consdata, y, delta);
                varmap[yname] = y;
                binary_fixed.insert(yname);
                redundant_var.insert(yname);
			}
        }
        std::cout<<"the number of infeasible constraints is  "<<redundant_var.size()<<std::endl;
        if (redundant_var.size() == 0)
        {
            std::unordered_map<std::string, double> tmpmask = { {"mask",1.0} };
            return { tmpmask,tmpsolution };
        }

        SCIP_CONS* lb_cut;
        SCIPcreateConsBasicLinear(scip,&lb_cut, "lb_cut", 0, NULL, NULL, neighbor-binary_fixed.size()+1, SCIPinfinity(scip));
        for (auto fixed : binary_free)
        {
        	SCIPaddCoefLinear(scip, lb_cut, varmap.at(fixed), 1);
        }
        for (auto fixed : binary_fixed)
        {
			SCIPaddCoefLinear(scip, lb_cut, varmap.at(fixed), -1);
		}	        
        SCIPaddCons(scip, lb_cut);
        SCIPwriteOrigProblem(scip, "test.lp", NULL, NULL);
        SCIPsolve(scip);
        std::cout << "first solve status : " << SCIPgetStatus(scip) << " first solve obj : " << SCIPgetPrimalbound(scip) << std::endl;
        SCIP_SOL* endsol = SCIPgetBestSol(scip);
        std::unordered_map<std::string, double> solution;
        if (endsol != NULL) // 检查是否存在解
        {
            // 获取模型的所有变量
            SCIP_VAR** vars = SCIPgetVars(scip);
            int nvars = SCIPgetNVars(scip);



            // 遍历每个变量，获取其在解中的值
            for (int i = 0; i < nvars; ++i)
            {
                std::string var_name = SCIPvarGetName(vars[i]);
                var_name = var_name.substr(2);
                if (redundant_var.find(var_name) == redundant_var.end())
                {
                    //std::cout<< var_name << std::endl;
                    double var_value = SCIPgetSolVal(scip, endsol, vars[i]);
                    solution[var_name] = var_value;
                }
            }
        }
        else
        {
            std::cout << "No feasible solution found." << std::endl;
        }
        std::cout << "first solve size " << solution.size() << std::endl;

        SCIPchgLhsLinear(scip,lb_cut, -SCIPinfinity(scip));
        SCIPchgRhsLinear(scip, lb_cut, neighbor - binary_fixed.size());
		SCIPsolve(scip);
        std::cout<<"second solve status : "<<SCIPgetStatus(scip)<<" second solve obj : "<<SCIPgetPrimalbound(scip)<<std::endl;
        SCIP_SOL* endsol_copy = SCIPgetBestSol(scip);
        std::unordered_map<std::string, double> solution_copy;
        if (endsol_copy != NULL) // 检查是否存在解
        {
            // 获取模型的所有变量
            SCIP_VAR** vars = SCIPgetVars(scip);
            int nvars = SCIPgetNVars(scip);
            // 遍历每个变量，获取其在解中的值
            for (int i = 0; i < nvars; ++i)
            {
                std::string var_name = SCIPvarGetName(vars[i]);
                if (redundant_var.find(var_name) == redundant_var.end())
                {
                    double var_value = SCIPgetSolVal(scip, endsol_copy, vars[i]);
                    solution_copy[var_name] = var_value;
                }
            }
        }
        else
        {
            std::cout << "No feasible solution found." << std::endl;
        }
        std::cout << "second solve size " << solution_copy.size() << std::endl;
        SCIPwriteOrigProblem(scip, "test2.lp", NULL, NULL);
        SCIPfree(&scip);
        return { solution,solution_copy };
	}

}
