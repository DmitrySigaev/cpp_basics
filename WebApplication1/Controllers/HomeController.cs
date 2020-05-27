using Microsoft.AspNetCore.Mvc;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace WebApplication1.Controllers
{
    [Route("")]
    public class HomeController : Controller
    {
        private readonly GuidProvider _guidProvider;
        public HomeController(GuidProvider provider)
        {
            _guidProvider = provider;
        }

        [HttpGet]
        public IActionResult GetId() => Ok(_guidProvider.Id);
    }
}
