using System;

namespace WebApplication1
{
    public class GuidProvider
    {

        public Guid Id { get; }
        public GuidProvider()
        {
            Id = new Guid();
        }
    }
}